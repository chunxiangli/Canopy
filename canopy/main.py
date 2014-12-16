import sys
from canopy import msg_exit
if sys.version_info < (2,7) :
    msg_exit("Sorry: requires python version 2.7 or greater (but not python 3.x)")

import os, time, random, time, shutil
from canopy.config import Config, ConfigAndOptionParser, get_number_of_cpus
from canopy.logger import MESSENGER, get_logger
from canopy.tools import Alignment, MultiAlignments, format_alignment_result_in_directory
from canopy.utils import ALIGNER_CLASS, MERGER_CLASS, TREE_ESTIMATOR_CLASS, TRANSLATOR_CLASS, translate_data, write_xml
from canopy.job import jobQueue, mainWorker
from canopy.file_manage import TempFileManager, makeArchive, remove_files_from_dir, remove_files
from canopy.coestimate import CoEstimator
from canopy.tree import PhylogeneticTree, write_tree_to_file, robinson_foulds_distance, symmetric_difference

_LOG = get_logger(__name__)

_RunningJob = None	
def initial_configuration():
	user_config = ConfigAndOptionParser()
	user_config.read_from_commandline()

	if not user_config.get("main", False):
		MESSENGER.send_error('There is no parameter. For the usage of the tool, please run "canopy --help" or "canopy -h"!!!')
		sys.exit()

	dir_name = user_config["main"].get("work_directory", os.curdir)
	seq_path = os.path.realpath(os.path.expanduser(user_config["main"]["seq"]))
	user_config["main"]["seq"] = seq_path
	if not os.path.exists(seq_path):
		MESSENGER.send_error("Input sequences file doesn't exist or not a valid file.")
		raise ValueError("Input sequences file donesn't exist or not a valid file")

	input_dir = os.path.dirname(seq_path)
	job_name = user_config["main"].get("name", "job")
	if dir_name is None:
		Config.work_directory = os.path.join(input_dir, job_name)
	else:
		#if running on cluster, get the local temporary directory, then copy compressed files to archive directory
                if "TMPDIR" in os.environ:
                        dir_name = os.path.join(os.environ["TMPDIR"], dir_name)

                dir_name = os.path.realpath(dir_name)
		parent_dir, child_dir = os.path.split(dir_name)

		if not os.path.exists(parent_dir):
			raise OSError("The parent directory of the specified main work direcotry %s doesn't exist."%dir_name)
                elif not os.path.exists(dir_name):
                        os.makedirs(dir_name)

                Config.work_directory = os.path.join(dir_name, job_name)

	if os.path.exists(Config.work_directory):
		wdir = Config.work_directory	
		suffix = ''
		import random, string
		while os.path.exists(wdir+suffix):
			suffix = ''.join([random.choice(string.letters) for i in range(5)])
		Config.work_directory = wdir + suffix
		os.makedirs(Config.work_directory)
		MESSENGER.send_info("Specified work directory '%s' existed. Results will be stored in '%s'"%(wdir, Config.work_directory))
	else:
		os.makedirs(Config.work_directory)

	Config.log_stream = open(os.path.join(Config.work_directory, "%s.log"%job_name), 'w')

	MESSENGER.out_log_streams.append(Config.log_stream)
	MESSENGER.err_log_streams.append(Config.log_stream)
	if not "num_cpus" in user_config["main"]:
		user_config["main"]["num_cpus"] = 1 
	else:
		num_cpus = user_config["main"]["num_cpus"]
		if num_cpus > 1:
			real_num_cpus = get_number_of_cpus()
			if num_cpus > real_num_cpus:
				user_config["main"]["num_cpus"] = max([1, real_num_cpus-1])
				MESSENGER.send_warning("Only %d cpus are detected."%real_num_cpus)

	#save the configuration
	user_config["main"]["work_directory"] = Config.work_directory
	export_config_path = user_config["main"].get("config_export_path", "config.back") 
	export_config_path = user_config.write_to_file(export_config_path)
	MESSENGER.send_info("Configurations are saved as %s"%export_config_path)

	Config.main = user_config["main"]
	del user_config["main"]
	Config.user_config = user_config

	check_prank_option()

def check_prank_option():
	if Config.main.get("aligner", "prank") == "prank" or Config.main.get("merger", "prank") == "prank":
		prank_args = Config.user_config["prank"]("prank.args", "")
		_LOG.debug("filter prank args:%s"%prank_args)
		if not "-showtree" in prank_args:
			Config.user_config["prank"]["prank.args"] = prank_args + " -showtree"	

		_LOG.debug("prank args:%s"%Config.user_config["prank"]["prank.args"])

		Config.main["output_options"] = []

		if Config.main.get("showall", False):
			Config.main["output_options"] = ["-showanc", "-showevents", "-showxml"]
		else:
			if Config.main.get("showanc", False):
				Config.main["output_options"].append("-showanc")

			if Config.main.get("showevents", False):
				Config.main["output_options"].append("-showevents")

			if Config.main.get("showxml", False):
				Config.main["output_options"].append("-showxml")

def post_file_process():
	data_name = os.path.basename(Config.work_directory)
	target_dir = os.path.realpath(os.path.expanduser(Config.main["archive_directory"]))

	parent, dir_name = os.path.split(target_dir)
	if "" == dir_name:
		parent, dir_name = os.path.split(parent)

	if not os.path.exists(parent):
		raise OSError("The parent directory of the specified archive directory %s doesn't exist."%target_dir)
	elif not os.path.exists(target_dir):
                os.makedirs(target_dir)

	target = os.path.join(target_dir, data_name)
	if not os.path.exists(target_dir):
		os.makedirs(target_dir)

	import random, string
	while os.path.exists(target+".tgz"):
		target = os.path.join( target_dir, data_name + ''.join([random.choice(string.letters) for i in range(5)]))
	target = target+".tgz"
	MESSENGER.send_info("Results will be stored as '%s'"%(target))

	archive_path = makeArchive(Config.work_directory)	
	shutil.copy2(archive_path, target)
	shutil.rmtree(Config.work_directory)
	os.remove(archive_path)
		
def store_result(alignment, prefix, tree_str=None, score=None, name_map=None):
	file_path = None
	suffix = ".tree"

	if isinstance(alignment, MultiAlignments):
		file_path = []
		for name in alignment.names:
			file_name = os.path.join(Config.work_directory, "%s_%s.fas"%(prefix, name))
			file_path.append(file_name)
	else:
		file_name = os.path.join(Config.work_directory, "%s.fas"%(prefix))
		file_path = file_name

	alignment.write_to_path(file_path, name_map=name_map)

	def write_tree(t_str, mid=""):
		st_file = os.path.join(Config.work_directory, prefix+mid+suffix)
		_LOG.debug(st_file)
                if name_map:
                        phy_tree = PhylogeneticTree.read_from_string(t_str)
                        phy_tree.rename_leaf_names(name_map)
                        write_tree_to_file(phy_tree, st_file, schema=tree_schema)
                else:
                        write_tree_to_file(t_str, st_file, schema=tree_schema)

        if tree_str is not None:
                tree_schema = Config.main.get("tree_format", "newick")
                if tree_schema == "nexus":
                        suffix = ".nex"
		if isinstance(tree_str, dict):
			[write_tree(tree_str[name], '_'+name) for name in tree_str]
		else:
			write_tree(tree_str)

	def write_score(sc, mid=""):
		with open(os.path.join(Config.work_directory, "%s%s_score.txt"%(prefix, mid)),'w') as score_file:
                        score_file.write("%.10f"%sc)
                
	if score is not None:
		if isinstance(score, dict):
			[ write_score(score[name], '_'+name) for name in score ]
		else:
			write_score(score)

	return file_path

def create_tools(tmpFileM, need_translator):	
	def initial_tool(tool_name, tool_type, args=" "): 
		class_tool_list = [ALIGNER_CLASS, MERGER_CLASS, TREE_ESTIMATOR_CLASS, TRANSLATOR_CLASS] 
		tool_class = None
		try:
			tool_class = [tool for tool in class_tool_list[tool_type] if tool_name == tool.name][0]
		except IndexError, e:
			MESSENGER.send_error("%s is not recognized. Please choose tool from: %s."%(tool_name, " ".join([tool.name for tool in class_tool_list[tool_type]])))
			raise e

		tool_args = args
		option_name = "%s.args"%(tool_class.name)

		if tool_class.name == 'prank' and Config.main.get("align_datatype", None) == 'CODON':
			tool_args += " -codon "

		if option_name in Config.user_config[tool_class.name]:
			tool_args += Config.user_config[tool_class.name][option_name]

		if tool_class.name == 'mafft' and Config.user_config['mafft'].get('nolocal', False):
                        return tool_class(tmpFileM, cmd=Config.user_config[tool_class.name]["%s.path"%(tool_class.name)], args=tool_args, nolocal=True)
                else:   
                        return tool_class(tmpFileM, cmd=Config.user_config[tool_class.name]["%s.path"%(tool_class.name)], args=tool_args)
	
	initial_aligner_name = Config.main.get("initial_aligner", 'mafft')
	Config.initial_aligner = initial_tool(initial_aligner_name, 0)

	MESSENGER.send_info("Initial aligner:%s created successfully."%initial_aligner_name)

	aligner_name = Config.main.get("aligner", 'prank')
	if initial_aligner_name == aligner_name:
		Config.aligner = Config.initial_aligner
	
	else:
		Config.aligner = initial_tool(aligner_name, 0) 
	MESSENGER.send_info("Iterative aligner:%s created successfully."%aligner_name)

	merger_name = Config.main.get("merger", 'prank')
	Config.merger = initial_tool(merger_name, 1)
	MESSENGER.send_info("Merge tool:%s created successfully."%merger_name)

	tree_estimator_name = Config.main.get("tree_estimator", 'raxml')
	Config.tree_estimator = initial_tool(tree_estimator_name, 2) 
	MESSENGER.send_info("Tree estimate tool:%s created successfully."%tree_estimator_name)

	if tree_estimator_name != 'raxml':
		Config.main["replicate_num"] = 0

	if need_translator:
		translator_name = Config.main.get("translator", 'prank')
		if translator_name == initial_aligner_name:
			Config.translator = Config.initial_aligner
		elif translator_name == aligner_name:
			Config.translator = Config.aligner
		else:
			Config.translator = initial_tool(translator_name, 3)

def read_input_files(align_datatype):
	input_seqs = None
	is_multi_alignments = False
	initial_score = None
	initial_tree = None

	try:
		
		data_datatype = Config.main.get("datatype", 'DNA')
		input_format = Config.main.get("format", 'fasta')

		
		#reading sequences
		seq_path = Config.main["seq"]
		err_msg = "Input sequence file is empty."
		if os.path.isdir(seq_path):
			is_multi_alignments = True
			seq_path = [ "%s/%s"%(seq_path, f) for f in os.listdir(seq_path) if not f.startswith(".")]
			input_seqs = MultiAlignments(align_datatype)	
		elif os.path.isfile(seq_path):
			input_seqs = Alignment(align_datatype)

		name_map = input_seqs.read_from_path(seq_path, input_format, data_datatype)

		if input_seqs.is_empty():
			err_msg +="Please make sure the files' format are correctly specified."
			MESSENGER.send_error(err_msg)	
			raise ValueError(err_msg)

		if name_map.keys() == name_map.values():
			name_map = None	
		
		if not is_multi_alignments:	
			initial_phy_tree = None
			if "tree" in Config.main:
				initial_phy_tree = PhylogeneticTree.read_from_path(os.path.expanduser(Config.main["tree"]))
			elif input_format == "nexus":
				initial_phy_tree = PhylogeneticTree.read_from_path(seq_path, "nexus")

			if "score" in Config.main:
				initial_score = float(Config.main["score"])


			if initial_phy_tree is not None:
				#Also need to encode leaf names
				if name_map is not None:
					initial_phy_tree.rename_leaf_names(name_map, restore=False)
				initial_tree = initial_phy_tree.as_newick_string()
				#tell rooted or unrooted tree
				initial_phy_tree.resolve_polytomies()
				istring = initial_phy_tree.as_newick_string()
				if istring.count("(") == istring.count(","):
					Config.main["rooted"] = True
		#TODO:whether need to accept initial tree and initial score for multiple alignments

		return input_seqs, initial_tree, initial_score, is_multi_alignments, name_map

	except Exception as e:
		raise e


def main():
	start = time.time()
	global _RunningJob

	initial_configuration()

	tempFileManager = TempFileManager()
	_top_level_temp = tempFileManager.create_top_level_temp(parent=Config.work_directory)
	
	save_option = Config.main.get("save_option", "simple")
	delete_temps = True
        if "all" == save_option:
                delete_temps = False
	_LOG.debug("delete_temps:%s"%delete_temps)

	try:
		data_datatype = Config.main.get("datatype", "DNA")
		align_datatype = Config.main.get("align_datatype", None)
		need_translator = False
		if align_datatype is not None and "DNA" == data_datatype and align_datatype != data_datatype:
			need_translator = True

		create_tools(tempFileManager, need_translator)
		
		tree_model = Config.user_config[Config.tree_estimator.name].get("%s.model"%Config.tree_estimator.name, "")

		input_seqs, initial_tree, initial_score, is_multi_alignments, name_map = read_input_files(align_datatype)

		_LOG.debug("Read input:%.2gs"%(time.time()-start))

	        #Original input file
		saved_input_path = store_result(input_seqs, "input", initial_tree, name_map=name_map)
		input_transalte_path = None	
		if need_translator and align_datatype == "PROTEIN":
			input_translate_path = translate_data(Config.translator, saved_input_path, Config.work_directory)
			input_seqs.read_from_path(input_translate_path, data_type="PROTEIN", prefix="input_")
			MESSENGER.send_info("The input sequences have been translated into %s."%(align_datatype.lower()))

		num_cpus = Config.main["num_cpus"]
		mainWorker.set_max_num_workers(num_cpus) 
		mainWorker.start_workers(num_cpus)	

		#initial alignment
		init_start = time.time()
		if initial_tree is None:
			initial_temp = tempFileManager.create_subdir("guide_tree_estimation")
			MESSENGER.send_info("Initial alignment start...")
			initial_input_seqs = type(input_seqs)()
			if need_translator and align_datatype == "CODON":
				translate_path = translate_data(Config.translator, saved_input_path, Config.work_directory)
				initial_input_seqs.read_from_path(translate_path, data_type="PROTEIN", prefix="input_")
				if is_multi_alignments:
					remove_files(translate_path)
				else:
					os.remove(translate_path)
			else:
				initial_input_seqs = input_seqs

			if is_multi_alignments:
				initial_score = {}
				initial_tree = {}
				_RunningJob = []
				for name in initial_input_seqs.names:
					job  = Config.initial_aligner.create_job(initial_input_seqs[name],
										 id=name,
										 tmp_dir=initial_temp,
										 delete_temps=delete_temps)
					jobQueue.put(job, initial_input_seqs[name].get_num_taxa())
					_RunningJob.append(job)

				for job in _RunningJob:
					initial_input_seqs[job.id].update(job.get_result())

				_RunningJob = None
				MESSENGER.send_info("Initial alignment done.")

				MESSENGER.send_info("Initial tree estimation start...")
				for name in initial_input_seqs:
					tree_estimate_job = Config.tree_estimator.create_job(initial_input_seqs[name],
									id=name,
									tmp_dir=initial_temp,
									num_cpus=num_cpus,
									delete_temps=delete_temps,
									model=tree_model)
					_RunningJob = tree_estimate_job
					jobQueue.put(tree_estimate_job, initial_input_seqs[name].get_num_taxa())
					initial_score[name], initial_tree[name] = tree_estimate_job.get_result()
					_RunningJob = None	
				MESSENGER.send_info("Initial tree estimation done.")
			else:
				job  = Config.initial_aligner.create_job(initial_input_seqs,
									 tmp_dir=initial_temp,
									 delete_temps=delete_temps,
									 description="whole")
				_RunningJob = job
				jobQueue.put(job, input_seqs.get_num_taxa())
				initial_input_seqs.update(job.get_result())
				_RunningJob = None
				MESSENGER.send_info("Initial alignment done.")

				MESSENGER.send_info("Initial tree estimation start...")
				tree_start = time.time()
				job = Config.tree_estimator.create_job(initial_input_seqs,
							       tmp_dir=initial_temp,
							       num_cpus=num_cpus,
							       delete_temps=delete_temps,
							       model=tree_model)
				_RunningJob = job
				jobQueue.put(job, initial_input_seqs.get_num_taxa())
				initial_score, initial_tree = job.get_result()
				_RunningJob = None
				MESSENGER.send_info("Initial tree estimation done.")
			
			saved_initial_result_path = store_result(initial_input_seqs, "initial", initial_tree, initial_score, name_map)

			_LOG.debug("Initial done.%.2gs"%(time.time()-init_start))

			if not Config.main.get("test", False):
				initial_input_seqs = None

			if need_translator:
				back_translated_path = translate_data(Config.translator,
								      saved_initial_result_path,
								      Config.work_directory,
								      dna_path=saved_input_path)

		MESSENGER.send_info("Iterative coestimation start...")

		start_co = time.time()

		def coestimate_single(align, tree, score, name=""):		
			jobQueue.whole_size = align.get_num_taxa()
			global _RunningJob
			co_temp = tempFileManager.create_subdir("divide_and_merge%s"%('_'+name))

			coestimator = CoEstimator(Config.aligner,
						  Config.merger,
						  Config.tree_estimator,
						  align,
						  tree,
						  file_manager=tempFileManager,
						  score=score,
						  tmp_dir=co_temp,
						  model=tree_model,
						  delete_temps=delete_temps,
						  translator=getattr(Config, "translator", None),
						  num_cpus=num_cpus,
						  max_prob_size=Config.main.get("max_prob_size", align.get_num_taxa()),
						  need_sub_iter=Config.main.get("need_sub_iter", False),
						  save_option=save_option, 
						  without_guide=Config.main.get("without_guide",False),
						  decomposition=Config.main.get("decomposite_strategy","seed"),
						  align_datatype=align_datatype,
						  tree_format=Config.main.get("tree_format", "newick"),
						  output_options=Config.main.get("output_options",[]),
						  id=name,
						  name_map=name_map,
						  rooted=Config.main.get("rooted", False),
						  replicate_num=Config.main.get("replicate_num", 0))

			_RunningJob = coestimator
			coestimator.start()	
			result_path, result_tree, result_score = coestimator.get_result()
			_RunningJob = None

			return result_path, result_tree, result_score

		writeXML = lambda: "showxml" in Config.main or "showall" in Config.main

		saved_result_path = None
		best_tree = None
		best_score = None
		if is_multi_alignments:
			saved_result_path={}
			best_tree = {}
			best_score = {}
			for name in input_seqs:
				saved_result_path[name], best_tree[name], best_score[name] = coestimate_single(input_seqs[name], initial_tree[name], initial_score[name], name)	
			#infer species tree
			input_seqs.read_from_path(saved_result_path.values(), data_type=input_seqs.datatype)
			MESSENGER.send_info("Species tree inference start...")
			concatenate_result = input_seqs.concatenate()
			species_temp = tempFileManager.create_subdir("species_tree_estimator")
			tree_estimate_job = Config.tree_estimator.create_job(concatenate_result,
							tmp_dir=species_temp,
							num_cpus=num_cpus,
							delete_temps=delete_temps,
							model=tree_model)
			_RunningJob = tree_estimate_job
			jobQueue.whole_size = None
			jobQueue.put(tree_estimate_job, concatenate_result.get_num_taxa())
			best_score, best_tree = tree_estimate_job.get_result()
			saved_concatenate_result = store_result(concatenate_result, "species", best_tree, best_score, name_map=name_map)

			if writeXML():
				write_xml(Config.work_directory, "species", concatenate_result, best_tree, name_map=name_map)

			concatenate_result = None
			_RunningJob = None 
                        MESSENGER.send_info("Species tree inference done.")
	
		else:
			saved_result_path, best_tree, best_score = coestimate_single(input_seqs, initial_tree, initial_score)
			if writeXML():
				input_seqs.read_from_path(saved_result_path)
				write_xml(Config.work_directory, "result_final", input_seqs, best_tree, name_map=name_map)

                #If the alignment target is PROTEIN while the input sequences are in DNA, we need back translate the final result into DNA.
		translate_result_path = None
                if need_translator and "PROTEIN" == align_datatype:
			_LOG.debug("need back translate")
			translate_file_name = "result_translated"
			if is_multi_alignments:
				saved_result_path = saved_result_path.values()	
				translate_file_name = "species_translated"

                        translate_result_path = translate_data(Config.translator, saved_result_path, Config.work_directory, dna_path=saved_input_path)

			input_seqs.read_from_path(translate_result_path, data_type="DNA", prefix="result_") 				
			input_seqs.write_to_path(os.path.join(Config.work_directory, "%s.fas"%translate_file_name), name_map=name_map)

			if writeXML():
				write_xml(Config.work_directory, translate_file_name, input_seqs, best_tree, name_map=name_map)

			_LOG.debug("finish back translate")

		if "output_format" in Config.main:
			_LOG.debug("need format alignment")
			filter_files = [saved_input_path]
			if is_multi_alignments:
				filter_files = saved_input_path

			if need_translator and align_datatype == "PROTEIN":
				if is_multi_alignments:
					filter_files.extend(translate_result_path)
				else:
					filter_files.append(translate_result_path)
				format_alignment_result_in_directory(Config.work_directory, Config.main.get("output_format"), "PROTEIN", filter_files)
			else:
				format_alignment_result_in_directory(Config.work_directory, Config.main.get("output_format"), data_datatype, filter_files)

			_LOG.debug("finish format alignment")

		canopy_span = time.time()-start_co

		_LOG.debug("Iterative coestimation Finished. Time:%.2gs"%canopy_span)

		if Config.main.get("test", False):
			two_phase(initial_tree, tempFileManager, input_seqs, delete_temps)


		#deal with tempfiles
		_LOG.info("save_option %s"%save_option)
                if 'all' == save_option:
                        tempFileManager.makeArchive()
		else:
                        tempFileManager.clear()

		if 'simple' == save_option:
			remove_files_from_dir(Config.work_directory, "iteration")
	except KeyboardInterrupt:
		MESSENGER.send_error("KeyboardInterrupt.")
	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		MESSENGER.send_error("At main stage. Msg:%s. Type:%s, FileName:%s, Line:%d"%(str(e), exc_type, fname, exc_tb.tb_lineno))
	finally:
		if "archive_directory" in Config.main:	
			MESSENGER.send_info("post_file")
			post_file_process()

		total_time = time.time() - start
		MESSENGER.send_info("finished. Total time:%.2gs"%total_time)
		try:
			Config.log_stream.close()
		except:
			pass
	
def two_phase(initial_tree, tempFileManager, alignment, delete_temps):
	if initial_tree is not None and not Config.main.get("rooted", False):
                phy_tree = PhylogeneticTree.read_from_string(initial_tree)
                phy_tree.reroot_at_midpoint(update_splits=True)
                phy_tree.resolve_polytomies()
                initial_tree = phy_tree.as_newick_string()
	
	MESSENGER.send_info("Two phase sequence alignment and phylogenetic tree construction start..")
	num_cpus = Config.main["num_cpus"]

	prank_args = " "
        if "prank.args" in Config.user_config['prank']:
                prank_args = Config.user_config['prank']['prank.args']
	from utils import PRANK
        prank = Prank(tempFileManager, cmd=Config.user_config['prank']['prank.path'], args=prank_args)
	tree_model = Config.user_config[Config.tree_estimator.name].get("%s.model"%Config.tree_estimator.name, "")

	prank_start = time.time()
	prank_temp = tempFileManager.create_subdir("twophase_prankGT")
	if Config.main.get("without_guide", False):
		prank_alignment = prank.run(alignment,
					    description="prank_align_without_guide_tree",
					    iterate=Config.main.get("max_iter",1),
                                            tmp_dir=prank_temp,
                                            delete_temps=False)
	else:
		prank_alignment = prank.run(alignment,
					    description="prank_align_with_guide_tree",
					    iterate=Config.main.get("max_iter",1),
                                            guide_tree=initial_tree,
                                            tmp_dir=prank_temp,
                                            delete_temps=False)

        log_file = os.path.join(prank_temp, "_temp_prank", "stdout.txt")
	last_score = None
	with open(log_file) as lfo:
        	results = filter(lambda x: "Alignment score:" in x, lfo)
        	last_score = results[-1].split(':')[1].strip()

	prank_score, prank_tree = Config.tree_estimator.run(prank_alignment, 
							    description="tree_estimate_on_prank_alignment",
                                                            tmp_dir=prank_temp,
                                                            num_cpus=num_cpus,
							    model=tree_model,
                                                            delete_temps=delete_temps)  

	store_result(prank_alignment,  'twophase_prankGT', prank_tree)
	if Config.tree_estimator.name != "pranktree":
                        with open(os.path.join(Config.work_directory, "twophase_prankGT_prank_score.txt"), 'w') as pscore:
                                pscore.write("%.5f %s"%(prank_score, last_score))

	prank_span = time.time() - prank_start  

	MESSENGER.send_info("PRANKGT time:%.2g"%prank_span)

	clustalw_args = ""
        if 'clustalw.args' in Config.user_config['clustalw']:
                clustalw_args = Config.user_config['clustalw']['clustalw.args']

	from utils import ClustaW
        clustalw = ClustalW(tempFileManager, cmd=Config.user_config['clustalw']['clustalw.path'])

	clustalw_start = time.time()
	clustalw_temp = tempFileManager.create_subdir("twophase_clustalw")
	clustalw_alignment = clustalw.run(alignment, description="clustalw_align", tmp_dir=clustalw_temp, delete_temps=delete_temps)

	clustalw_score, clustalw_tree = Config.tree_estimator.run(clustalw_alignment,
                                                                  description="tree_estimate_on_clustalw_alignment",
                                                                  tmp_dir=clustalw_temp,
                                                                  num_cpus=num_cpus,
								  model=tree_model,
                                                                  delete_temps=delete_temps)   

	store_result(clustalw_alignment,  'twophase_clustalw', clustalw_tree)
	clustalw_span = time.time() - clustalw_start
	MESSENGER.send_info("Clustalw time:%.2g"%clustalw_span)

if __name__ == "__main__":
	main()
