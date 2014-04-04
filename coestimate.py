import os, time, sys, shutil
from StringIO import StringIO
from threading import Event
from dendropy import Tree
from job import jobQueue, JobBase
from tree import PhylogeneticTree, write_tree_to_file, reroot_at_midpoint, robinson_foulds_distance, symmetric_difference
from logger import get_logger, MESSENGER
from config import Config
from tools import MultiAlignments, Alignment, alignment_accuracy
from file_manage import copy_files, remove_files
from utils import translate_data,generate_prank_output

_LOG = get_logger(__name__)
_DEFAULT_MAX_ITER = 5

class CoEstimator(JobBase):
	def __init__(self, aligner, merger, tree_estimator, alignment, tree, **kwargs):
		JobBase.__init__(self, **kwargs)
		self._alignment = alignment
                self._num_taxa = alignment.get_num_taxa()
		self._tree = tree
		self._aligner = aligner
		self._merger = merger
                self._align_datatype = kwargs.get("align_datatype", None)
		self._tree_estimator = tree_estimator 
		self._translator = kwargs.get("translator", None) 
		self._file_manager = kwargs.get("file_manager")
		self._save_option = kwargs.get("save_option", "simple")
		self._best_score = kwargs.get("score", None)
		self._subiter = False
		self._best_tree = self._tree
		self._best_iter = 0 
		self._new_tree = None
		self._old_tree = None
		self._event_list = Event()
		self._job = None
		self._killed = False
		self._iterational_result_files = []

	def _keep_iterating(self):	
		self.check_status()
		if 0 == self._cur_iter:
			self._start_time = time.time()
		if self._num_taxa < 4 and self._cur_iter > 1:
			return False

		max_time = int(Config.main.get("max_time", 0))
		if max_time > 0:
			if (time.time() - self._start_time) > max_time:
				return False
			else:
				return True
		elif self._cur_iter > int(Config.main.get("max_iter", _DEFAULT_MAX_ITER))-1:
			return False
		else:
			return True

	def need_split(self):
		return self._kwargs.get("max_prob_size", self._num_taxa) < self._num_taxa

	def start(self):	
		if self._kwargs.get("need_sub_iter") is None:
			self._subiter = True

		self.run_time += 1

		assert self._tree_estimator is not None

		self._workdir = self._kwargs.get("tmp_dir", os.curdir)
		self._cur_iter = 0
		self._start_time = time.time()
		self._num_cpus = self._kwargs.get("num_cpus", 1)

		self._name_encode = False
		if 0 == self._cur_iter:
			self._name_map = self._kwargs.get("name_map")
			if self._name_map:
				self._name_encode = True 

		k = dict(self._kwargs)
		k["score"] = None
		k["num_cpus"] = 2

		try:
			while self._keep_iterating():
				cur_iter_work_tmp_dir = self._file_manager.create_temp_subdir(parent=self._workdir, dir_name="step%d"%self._cur_iter)
				k["tmp_dir_parent"] = cur_iter_work_tmp_dir
				k["old_tree"] = self._old_tree

				MESSENGER.send_info("Iteration %d: align start..."%(self._cur_iter))

				self.align(**k)

				MESSENGER.send_info("Iteration %d: align done."%(self._cur_iter))
				
				#estimate new tree and check whether got improved
   				new_score = None
                                if self._num_taxa > 3:
                                        MESSENGER.send_info("Iteration %d: tree estimation start..."%(self._cur_iter))
                                        new_score = self._update_tree(cur_iter_work_tmp_dir)
                                        MESSENGER.send_info("Iteration %d: tree estimation done."%(self._cur_iter))
                                        if not self._subiter:
                                                self._store_iterational_tree_result(new_score)

                                        self._tree = self._new_tree
                                else:   
                                        _LOG.info("Can't use %s. The alignment size is smaller than 4."%self._tree_estimator.name)
                                        self._best_iter = self._cur_iter

				#self._tree = self._best_tree

				if "all" != self._save_option:
					MESSENGER.send_info("Iteration %d: remove current work directory."%self._cur_iter)
					self._file_manager.remove_dirs([cur_iter_work_tmp_dir])

				self._cur_iter += 1

			if not self._subiter:
				self._store_optimal_result()

		except Exception as e:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			raise RuntimeError("Unexpected error:%s.Type:%s, FileName:%s,Line:%d"%(str(e), exc_type, fname, exc_tb.tb_lineno))
		finally:
			self._event_list.set()

	def _update_tree(self, work_directory):
		new_concatenated_alignment = self._alignment
		new_alignment_num_taxa = self._alignment.get_num_taxa()

		assert new_alignment_num_taxa == self._num_taxa, "The number of taxa of new alignment not equals the original number."
		
		_LOG.debug("Start tree estimation...")
		self.check_status()
		
		if self._align_datatype == "CODON":
			_LOG.debug("Codon alignment need translation for tree estimation...")
			temp_result_file = "%s/iteration%d_temp_result.fas"%(self._workdir, self._cur_iter)
			new_concatenated_alignment.write_to_path(temp_result_file)

			new_concatenated_alignment = Alignment()
			translated_path = translate_data(self._translator, temp_result_file, self._workdir)
			new_concatenated_alignment.read_from_path(translated_path, data_type="PROTEIN")
			remove_files([temp_result_file, translated_path])
			_LOG.debug("translation done.")

		tree_estimator = self._tree_estimator.create_job(new_concatenated_alignment,
								 num_cpus = self._num_cpus,
								 tmp_dir=work_directory, 
								 delete_temps=self._kwargs.get("delete_temps"),
								 model=self._kwargs.get("model", ""),
								 id=self.id,
								 description=self._kwargs.get("description", "whole"))
		self._job = tree_estimator
		self.check_status()
		tree_estimator.start()
		new_score, self._new_tree = tree_estimator.get_result()
		new_concatenated_alignment = None
		self._job = None

		_LOG.debug("Tree estimation done.")
		def score_improved(new_score):

			if self._tree_estimator.name == "prankTree":
				if self._best_score is None:
					self._best_score = self._prank_score
					return True
				elif self._prank_score < self._best_score:
					self._best_score = self._prank_score
					return True
			elif self._best_score is None or new_score > self._best_score:
				self._best_score = new_score
				return True
			return False

		self.check_status()

		if score_improved(new_score):
			MESSENGER.send_info("Tree improved at iteration %d."%self._cur_iter)
			self._best_tree = self._new_tree
			self._best_iter = self._cur_iter
		
		#print "get new tree done"
                return new_score


	def align(self, **kwargs):
		phy_tree = PhylogeneticTree.read_from_string(self._tree)
        	phy_tree.reroot_at_midpoint(update_splits=True)
        	phy_tree.resolve_polytomies()

		MESSENGER.send_info("Create align job")
		self._create_align_job_for_single_data(self._alignment, phy_tree, **kwargs)

		self.check_status()

		_LOG.debug("Start get align result...")

		self._alignment.update(self._job.get_result())	

		#store iterational align results
		self._store_iterational_align_result()
		self._job = None
		self._old_tree = phy_tree.as_newick_string()

		_LOG.debug("Align done.")


	def _create_align_job_for_single_data(self, alignment, phy_tree, **kwargs):
                dir_parent = kwargs.get("tmp_dir_parent", os.curdir)
                max_prob_size = kwargs.get("max_prob_size", alignment.get_num_taxa())
                job = None

		dir_name = "align"
                kwargs["id"] = self.id

                align_work_dir = self._file_manager.create_temp_subdir(dir_parent, dir_name)
		kwargs["tmp_dir"] = align_work_dir
                if max_prob_size < phy_tree.num_taxa():
			_LOG.debug("Start divide and merge...")
			kwargs["datatype"] = alignment.datatype
			kwargs["description"] = kwargs.get("description", "whole")
			job = self._divide_and_merge(alignment, phy_tree, **kwargs)
			_LOG.debug("End divide and merge.")
                else:
                        guide_tree = phy_tree.as_newick_string()
			old_tree = self._old_tree
                        if kwargs.get("without_guide", False):
                                guide_tree = None
				old_tree = None

                        job = self._aligner.create_job(alignment, 
						       guide_tree=guide_tree,
                                                       old_tree=old_tree,
                                                       tmp_dir=align_work_dir,
                                                       delete_temps=False,
                                                       output_options=kwargs.get("output_options", []),
                                                       id=self.id,
                                                       description=kwargs.get("description", "whole"))
			if not self._subiter:
				jobQueue.put(job)
			else:
				try:
					job.start()
					job.get_result()
				except Exception as e:
					raise e

                self._job = job

	def _store_iterational_align_result(self):
		file_suffix = "fas"
		result_dir = Config.work_directory
		if self._subiter:
			result_dir = self._workdir

		output_options = self._kwargs.get("output_options",[])
		job = self._job
		wdir = job.cwd
		file_prefix = "iteration%d"%(self._cur_iter)
		if job.id != "":
			file_prefix = "iteration%d_%s"%(self._cur_iter, job.id)

		result_align_path = None

		if self._aligner.name == "prank" and self._merger.name.startswith("prank"):
			MESSENGER.send_info("need to deal with prank result")
			if self._name_encode:
				_LOG.debug("name encoded")
				result_align_path = generate_prank_output(result_dir,
									wdir,
									self._name_map,
									datatype=self._alignment.datatype,
									output_options=output_options,
									file_prefix=file_prefix)
			else:
				source_file_prefix = "result.aligned.best"
				if "".join(job._command).find("-ot=") != -1:
					source_file_prefix = "result.aligned"

				if -1 != job.name.find("prankmerger"):
					source_file_prefix = "merged"

				copy_files(wdir, source_file_prefix, file_prefix, result_dir)
				result_align_path = os.path.join(result_dir, file_prefix+".fas")

			if not self._subiter:	
				_LOG.debug("start to store prank score")
				log_file = os.path.join(wdir, "stdout.txt")
				results = filter(lambda x: -1 != x.find("Alignment score:"), open(log_file).readlines())
				self._prank_score = results[-1].split(':')[1].strip()
				_LOG.debug("prank score:%s"%self._prank_score)

			if 'all' != self._save_option:
				_LOG.debug("remove %s"%wdir)
				self._file_manager.remove_dirs([wdir])

			_LOG.debug("finish store prank score")
		else:
			_LOG.debug("Just output the alignment result")
			result_align_path = os.path.join(result_dir, "%s.%s"%(file_prefix, file_suffix))
			self._alignment.write_to_path(result_align_path, name_map=self._name_map)
			if self._best_score is not None:
				open(os.path.join(Config.work_directory, "%s_score.txt"%file_prefix), 'w').write("%.5f"%(self._best_score))


		self._iterational_result_files.append(result_align_path)

		_LOG.debug("result_aligns:%s"%result_align_path)

	def _store_optimal_result(self):
		MESSENGER.send_info("save final result...")

		work_dir = Config.work_directory		
		if self._subiter:
			work_dir = self._workdir

		source_prefix = "iteration%d"%self._best_iter
		optimal_prefix = "result"
		copy_files(work_dir, source_prefix, optimal_prefix)

	def _divide_and_merge(self, alignment, phy_tree, **kwargs):
                work_dir = kwargs.get("tmp_dir", os.curdir)
		#generate the AlignMergeTree
		_ALIGN_JOBS = []
	
		_LOG.debug("Start split...")
		alignMergeTree = AlignMergeTree(phy_tree, alignment, work_dir, self._aligner, self._merger, self._tree_estimator, **kwargs)
		_LOG.debug("End split.")

		if self._killed:
				[ job.kill() for job in _ALIGN_JOBS]
				raise RuntimeError("Coestimator Killed during align.")

		MESSENGER.send_info("The size of merge list:%d"%len(_MERGE_NODES))
		_merge_jobs = []
		_LOG.debug("Start merge...")
		while _MERGE_NODES:
			update = False
			if self._killed:
				[ job.kill() for job in _merge_jobs]
				raise RuntimeError("Coestimator Killed during merging.")

			for node in _MERGE_NODES:
				job = node._create_merge_job()
				if job:
					update = True
					jobQueue.put(job)
					_merge_jobs.append(job)
					_MERGE_NODES.remove(node)
					if not _MERGE_NODES:#return the root job
						return job
		_LOG.debug("End merge...")

	def _store_iterational_tree_result(self, new_score):
                _LOG.info("store_iterational_tree_result")
		suffix = ".tre"
		tree_schema = self._kwargs.get("tree_format", "newick")
		wdir = Config.work_directory
		if self._subiter:
			wdir = self._workdir

                if tree_schema == "nexus":
                        suffix = ".nex"

		score_file = os.path.join(wdir, "iteration%d_score.txt"%self._cur_iter)
		if self.id != "":
			score_file = os.path.join(wdir, "iteration%d_%s_score.txt"%(self._cur_iter, self.id))
			
                with open(score_file, 'w') as sf:
                        if new_score is not None:
                                sf.write("%.5f "%new_score)
                        try:   
                                sf.write("%s"%self._prank_score)
                        except:
                                pass

		if self._new_tree is not None:
			st_file = os.path.join(wdir, "iteration%d%s"%(self._cur_iter, suffix))
			if self.id != "":
				st_file = os.path.join(wdir, "iteration%d_%s%s"%(self._cur_iter, self.id, suffix))
				
			if self._name_encode:
				phy_tree = PhylogeneticTree.read_from_string(self._new_tree)
				phy_tree.rename_leaf_names(self._name_map)
				write_tree_to_file(phy_tree, st_file, schema=tree_schema)
			else:   
				write_tree_to_file(self._new_tree, st_file, schema=tree_schema)
	
	def wait(self):
		self._event_list.wait()

	def get_result(self):
		self.wait()
		if self._subiter:
			best_alignment = Alignment()
			best_alignment.read_from_path(self._iterational_result_files[self._best_iter], data_type=self._alignment.datatype)
			self._result = (best_alignment, self._best_tree, self._best_score)
		else:
			best_keyword = "iteration%d"%self._best_iter
			_LOG.debug("get_result:%s"%best_keyword)
			_LOG.debug("best result:%s"%self._iterational_result_files[self._best_iter])
			fname = self._iterational_result_files[self._best_iter]
			result_file = fname.replace(best_keyword, "result")
			self._result = (result_file, self._best_tree, self._best_score)
			MESSENGER.send_info("Best iter:%d"%self._best_iter)
		return self._result

	def kill(self):
		self._killed = True
		if self._job is not None:
			self._job.kill()

	def check_status(self):
		if self._killed:
			raise RuntimeError("Coestimator Killed.")
_MERGE_NODES = []
_ALIGN_JOBS = []
class AlignMergeTree(object):
        def __init__(self, phy_tree, alignment, work_dir, aligner, merger, tree_estimator, **kwargs):
                self._phy_tree = phy_tree
		self._num_taxa = self._phy_tree.num_taxa()
                self._tree = phy_tree.as_newick_string()
		self._alignment = alignment
		self.aligner = aligner
		self.merger = merger
		self.tree_estimator = tree_estimator
                self._id = kwargs.get("id","")
                self._align_job = None
                self._merge_job = None
                self._rChild = None
                self._lChild = None
                self._result = None
                self._kwargs = kwargs
                self._work_dir = work_dir
		self.generate_children()

	def get_result(self):
		if self._align_job is not None:
			self._align_job.check_status()
			self._result = self._align_job.result
		elif self._merge_job is not None: 
			self._merge_job.check_status()
			self._result = self._merge_job.result
		return self._result


        @property
        def id(self):
                return self._id

        @property
        def tree(self):
                return self._tree

	def need_split(self):
		return self._kwargs.get("max_prob_size", self._num_taxa) < self._num_taxa

	def generate_children(self):
		if self.need_split():
			if "seed" == self._kwargs.get("decomposition", "seed"):
				t1, t2, self.merge_dist1, self.merge_dist2 = self._phy_tree.bipartition_by_seed()
			else:
				t1, t2, self.merge_dist1, self.merge_dist2 = self._phy_tree.bipartition_by_longest_internal_edge()	
			
			k = dict(self._kwargs)
			file_manager = self._kwargs.get("file_manager")
			rChild_wdir = file_manager.create_temp_subdir(self._work_dir, '0')
			a1 = self._alignment.sub_alignment(t1.leaf_node_names())
			k["description"] += '/0'
			self._rChild = AlignMergeTree(t1, a1, rChild_wdir, self.aligner, self.merger, self.tree_estimator, **k)

			lChild_wdir = file_manager.create_temp_subdir(self._work_dir, '1')
			a2 = self._alignment.sub_alignment(t2.leaf_node_names())
			k["description"] += '/1'
			self._lChild = AlignMergeTree(t2, a2, lChild_wdir, self.aligner, self.merger, self.tree_estimator, **k)
			_MERGE_NODES.append(self)
		else:
			_LOG.debug("subalign create job")
			job = self._create_align_job()
			jobQueue.put(job)

	def _create_align_job(self):
		file_manager = self._kwargs.get("file_manager")
                job = None
                _LOG.info("Subalign job:%s creating..."%self._kwargs.get("description"))

                if self._num_taxa > 3 and self._kwargs.get("need_sub_iter", False):
                        k = dict(self._kwargs)
                        del k["need_sub_iter"]
                        del k["max_prob_size"]
			del k["output_options"]
			del k["name_map"]
                        k["tmp_dir"] = self._work_dir
			if "score" in k:
				del k["score"]

                        job = CoEstimator(self.aligner, self.merger, self.tree_estimator, self._alignment, self._tree, **k)
                else:
                        guide_tree = self._tree
			old_tree = self._kwargs.get("old_tree", None)
                        if self._kwargs.get("without_guide", False):
                                guide_tree = None
				old_tree = None

                        job = self.aligner.create_job(self._alignment,
			                              guide_tree,
                                                      old_tree=old_tree,
                                                      id=self.id,
                                                      tmp_dir=self._work_dir,
                                                      delete_temps=self._kwargs.get("delete_temps"),
                                                      description=self._kwargs.get("description"))
                self._align_job = job
               	_LOG.info("Subalign job:%s created..."%self._work_dir)

                return job

        def _create_merge_job(self):
                result1 = self._rChild.get_result()
                result2 = self._lChild.get_result()

		if result1 is None or result2 is None:
               		#_LOG.info("Children results are not ready.%s"%self._work_dir)
			return False
               	#_LOG.info("Children results are ready.%s"%self._work_dir)
                a1 = None
                a2 = None
                t1 = None
                t2 = None

                #When subalig job needs iteration, its result is a tuple
                if isinstance(result1, tuple):
                        a1, t1 = result1[:2]
			t1 = reroot_at_midpoint(t1)
                else:
                        a1 = result1
                        t1 = self._rChild.tree

                _LOG.debug("%s rchild result:%d"%(self._rChild._kwargs.get("description"), a1.get_num_taxa()))
		a1.dna_freqs = self._rChild._alignment.dna_freqs
                self._rChild = None

                if isinstance(result2, tuple):
                        a2, t2 = result2[:2]
			t2= reroot_at_midpoint(t2)
                else:
                        a2 = result2
                        t2 = self._lChild.tree
                _LOG.debug("%s lchild result:%d"%(self._lChild._kwargs.get("description"), a2.get_num_taxa()))
		a2.dna_freqs = self._lChild._alignment.dna_freqs
                self._lChild = None

		_LOG.debug("submerge:delete_temps:%s"%self._kwargs.get("delete_temps"))
		merge_tree = '(t1:{0},t2:{1})'.format(self.merge_dist1, self.merge_dist2)
		if self._kwargs.get("without_guide", False):
			t1 = None
			t2 = None
			merge_tree = None

		job = self.merger.create_job(a1, a2,
					     guide_tree1=t1,
                                             guide_tree2=t2,
		                             tmp_dir=self._work_dir,
                                             id=self.id,
                                             tree=self._tree,
                                             merge_tree=merge_tree,
                                             delete_temps=False,
                                             output_options=self._kwargs.get("output_options",[]),
                                             description=self._kwargs.get("description"))

                self._merge_job = job

               	_LOG.info("Merge job created...%s"%self._work_dir)
                return job

        def get_alignment(self, whole_alignment):
                return whole_alignment.sub_alignment(self._phy_tree.leaf_node_names())
