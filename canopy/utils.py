import sys, os, time, random, platform, shlex, re
from cStringIO import StringIO
from canopy.tree import PhylogeneticTree, write_tree_to_file
from canopy.logger import get_logger
from canopy.config import DEBUG
from canopy.job import Job, FakeJob, mainWorker, jobQueue
from canopy.tools import Alignment, MultiAlignments

_LOG = get_logger(__name__)

def file_checker(file_path):
	if not file_path:
		return False, "Path:%s is empty."%(file_path)

	if os.path.isfile(file_path):
		return True, file_path
	else:
		return False, "%s is not existed."%file_path

def translate_data(translator, seq_path, wdir, dna_path=None):
	assert translator is not None, "Please specify a tool to translate the data."

	args = ""
	backtranslate = False
	new_path = None

	if dna_path is not None:
		backtranslate = True
		
	if isinstance(seq_path, list):
		new_path = []
		for path in seq_path:
			if backtranslate:
				file_common_substr = "_".join(os.path.basename(path).split("_")[1:])
				dnafile = filter(lambda x: file_common_substr in x, dna_path)[0]
				args = "-dna=%s"%dnafile
			translated_file_path = translator.translate(path, wdir, args=args)
			new_path.append(translated_file_path)
	else:
		if backtranslate:
			args = "-dna=%s"%dna_path
		new_path = translator.translate(seq_path, wdir, args=args)	

	return new_path

def read_internal_alignment(result_file, schema="fasta", datatype=None, del_dirs=[], file_manager=None):
	aligned = Alignment()
	aligned.read_from_path(result_file, file_format=schema, data_type=datatype)

	if len(aligned) >= 1:
		#keep the intermediat file when error happend
		if len(del_dirs):
			file_manager.remove_dirs(del_dirs)
		return aligned
	else:
		raise ValueError("Program Failed:The alignment file %s is empty."%result_file)
	
def read_raxml_results(workdir, del_dirs, file_manager):
        file_list = os.listdir(workdir)
        tree_id = None
	log_score = None
	tree_str = ""
        for f in file_list:
                if f.startswith("RAxML_log"): 
                        tree_id = f.split(".")[1]
                        break
        raxml_log = os.path.join(workdir, "RAxML_log.%s"%tree_id)
	with open(raxml_log, "rU") as log_file:
        	log_score = filter(lambda x: "Final GAMMA_based Score of best tree" in x, log_file) 

        score = None
        if log_score:
                score = float(log_score[0].split(" ")[-1])
        else:   
		with open(raxml_log, "rU") as log_file:
                	score = float(log_file.readlines()[-1].split()[1])

        raxml_result = os.path.join(workdir, "RAxML_result.%s"%tree_id)
	with open(raxml_result, "rU") as tree_file:
        	tree_str = tree_file.read().strip()

	if len(del_dirs):
        	file_manager.remove_dirs(del_dirs)

        return score, tree_str

def read_fasttree_results(result_file, workdir, del_dirs, file_manager):
	tree_str = ""
	log_score = None
	with open(result_file, 'r') as tree_file:
		tree_str = re.sub("\)[0-9.]*:","):", tree_file.read(-1).strip("\n"))

	with open(workdir+"/stderr.txt", 'r') as log_file:
		log_score = filter(lambda x: "Optimize all lengths:" in x, log_file)

	score = None
	if log_score:
		score = float(log_score[0].split('=')[1].split()[0])
	if len(del_dirs):
		file_manager.remove_dirs(del_dirs)

	return score, tree_str
	

_internal_name_pattern = re.compile(r"\)[0-9\.]*:")
def read_phyml_results(workdir, del_dirs, file_manager):
        file_list = os.listdir(workdir)
        prefix_id = None
	log_score = None
	tree_str = ""

        for f in file_list:
		if "stats" in f:
                        prefix_id = f.split("stats")
                        break

	with open(os.path.join(workdir, "{0}stats{1}".format(prefix_id[0], prefix_id[1])), 'rU') as log_file:
		log_score = filter(lambda x: "Log-likelihood:" in x, log_file)

        score = None
        if log_score:
                score = float(log_score[0].split(":")[-1])

	with open(os.path.join(workdir, "%stree%s"%(prefix_id[0], prefix_id[1])), "rU") as tree_file:
        	tree_str = tree_file.read().strip()

	tree_str = re.sub(_internal_name_pattern, "):", tree_str)
	if len(del_dirs):
        	file_manager.remove_dirs(del_dirs)

        return score, tree_str

def generate_prank_output(target_dir, source_dir, name_map, datatype, output_options=[], output_format="fasta", file_prefix=None):
	_LOG.debug("start to generate PRANK output:%s"%output_options)
       	file_list = os.listdir(source_dir)
        prefix = None
	middle = ""
	restore = True
	if not isinstance(name_map, dict):
		restore=False

        for f in file_list:
		if "anc.best.fas" in f:
                        prefix = f.split("anc.best.fas")[0]
			middle = "best."
		elif "best.fas" in f:
                        prefix = f.split("best.fas")[0]
			middle = "best."
                elif "anc.fas" in f:
                        prefix = f.split("anc.fas")[0]
                        break
                elif f.endswith("fas"):
                        prefix = f.split("fas")[0]

        if file_prefix is None:
                file_prefix = prefix[:-1]

        #default
        origin_align_path = os.path.join(source_dir, prefix+middle+"fas")
        result_align_path = os.path.join(target_dir, "%s.fas"%file_prefix)
        out_anc = False

        if "-showanc" in output_options:
                out_anc = True
                origin_align_path = os.path.join(source_dir, prefix+"anc."+middle+"fas")

        alignment = Alignment()
        alignment.read_from_path(origin_align_path, file_format=output_format, data_type=datatype)
        alignment.write_to_path(result_align_path, file_format=output_format, name_map=name_map, suppress_ancester=out_anc)

        if out_anc:
                result_anc_align_path = os.path.join(target_dir, file_prefix+".anc.fas")
                alignment.write_to_path(result_anc_align_path, file_format=output_format, name_map=name_map)

        #tree
	result_tree_path = None
	result_anc_tree_path = None
        origin_tree_path = os.path.join(source_dir, prefix+middle+"dnd")
        result_tree_path = os.path.join(target_dir, file_prefix+".dnd")

	if out_anc:
                origin_tree_path = os.path.join(source_dir, prefix+"anc."+middle+"dnd")
                result_anc_tree_path = os.path.join(target_dir, file_prefix+".anc.dnd")
        ot = PhylogeneticTree.read_from_path(origin_tree_path)
	if restore:
		ot.rename_leaf_names(name_map)
	with open(result_tree_path, 'w') as rtp:
		ot.write(rtp, suppress_rooting=True, suppress_internal_node_labels=True)

	if out_anc:
	#anc tree
		with open(result_anc_tree_path, 'w') as ratp:
			ot.write(ratp, suppress_rooting=True)

        def replace_from_dict(str):
                if file_type == "events" and str.startswith("branch"):
                        prefix, name = str.strip().split()
                        if name in name_map:
                                str = '%s %s\n'%(prefix, name_map[name])
                elif file_type == "xml" and "name=" in str:
                                prefix, name = str.split('name="')
                                name = name.strip('">\n')
                                str = '%s name="%s">\n'%(prefix, name_map[name])
                return str
        if "-showevents" in output_options:
                #events
                origin_events_path = os.path.join(source_dir, prefix+middle+"events")
                result_events_path = os.path.join(target_dir, file_prefix+".events")

                file_type = "events"
                events_list = file(origin_events_path).readlines()
                if not out_anc:
			event_tree_str = events_list[4].strip()
			if event_tree_str == "":
				event_tree_str = events_list[3].strip()
                        ot = PhylogeneticTree.read_from_string(event_tree_str)
			if restore:
                        	ot.rename_leaf_names(name_map)
                raep_events = map(replace_from_dict, events_list[4:])
                with open(result_events_path, 'w') as rep:
                        rep.write("Alignment topology with node labels:\n\n%s\n"%ot.as_newick_string())
			rep.writelines(raep_events)
        if "-showxml" in output_options:
                #xml
                origin_xml_path = os.path.join(source_dir, prefix+middle+"xml")
                result_xml_path = os.path.join(target_dir, file_prefix + ".xml")
		_LOG.debug("result_xml_path:%s"%result_xml_path)
                file_type = "xml"
                with open(result_xml_path, 'w') as rxp:
                        rxp.writelines(map(replace_from_dict, file(origin_xml_path)))
	
	_LOG.debug("finish generating PRANK output")

	return result_align_path

def write_xml(directory, file_prefix, alignment, tree_str, name_map):
	if isinstance(alignment, MultiAlignments):
		alignment = alignment.concatenate()

	with open("%s/%s.xml"%(directory, file_prefix), 'w') as xmlFile:
		tree = PhylogeneticTree.read_from_string(tree_str, 'newick')
		xmlFile.write("<ms_alignment>\n<newick>\n")
		xmlFile.write(tree.as_newick_string())
		xmlFile.write("\n</newick>")
		xmlFile.write("<nodes>\n")
		for name in tree.leaf_node_names():
			if name_map is not None:
				xmlFile.write('<leaf id="%s" name="%s">\n'%(name, name_map[name]))	
			else:
				xmlFile.write('<leaf id="%s">'%name)	
			xmlFile.write(' <sequence>\n %s\n </sequence>\n</leaf>\n'%alignment[name])
		xmlFile.write("</nodes>\n</ms_alignment>\n")

class Tool(object):
	is_bundled = False

	def __init__(self, file_manager, **kwargs):
		self.file_manager = file_manager
		self.cmd = kwargs["cmd"]
		self.delete_temps = kwargs.get("delete_temps", True)

		self.user_config = None	
		args = kwargs.get("args")
		if args is not None and args != '':
			self.user_config = shlex.split(args)

		res, msg = self.check_executable()
		#if not os.path.exists(self.cmd):
		if not res:
            		if self.is_bundled:
                		err_msg = "The command '%s' does not exist. Please check the installation and try again." % self.cmd
            		else:
                		err_msg = "Executable file %s not found. Please install %s and/or configure its location correctly." % (self.cmd, self.name)
            		raise ValueError(err_msg)

	def check_executable(self):
		r, msg = file_checker(self.cmd)
		if not r:
			return (r, msg)
		else:
			if os.access(self.cmd, os.X_OK):
				return True, self.cmd
			else:
				return False, self.cmd

	def make_workdir(self, dir_name, prefix):
		assert dir_name is not None, "Must specify tmp_dir for %s."%self.name

		prefix = "%s_temp_%s"%(prefix, self.name)
		workdir = self.file_manager.create_temp_subdir(parent=dir_name, dir_name=prefix)

		return workdir

  	def create_job(self, *args, **kwargs):
        	raise NotImplementedError("Abstract Tool class cannot spawn jobs.")

	def run(self, *args, **kwargs):
		try:
			job = self.create_job(*args, **kwargs)
			job.start()
			return job.get_result()
		except Exception as e:
			raise e

class Aligner(Tool):
	result_prefix = "result.aligned"

	def __init__(self, file_manager, **kwargs):
		Tool.__init__(self, file_manager, **kwargs)	
		
	def _preprocess_input(self, alignment, **kwargs):
		workdir_path = self.make_workdir(kwargs.get("tmp_dir"), kwargs.get("id", ""))
        	input_file = os.path.join(workdir_path, "input.fasta")
        	alignment.write_to_path(input_file)
        	result_file = os.path.join(workdir_path, self.result_prefix)

        	return workdir_path, input_file, result_file

	def check_and_create_fake_job(self, alignment, **kwargs):
		num_taxa = alignment.get_num_taxa()
		if num_taxa < 2:
			kwargs["description"] = "%s_fake%s"%(kwargs.get("description", "whole"), self.name)
                        return True, FakeJob(alignment, **kwargs)

		return False, None
	

	
	def create_job(self, *args, **kwargs):
		raise NotImplementedError("Abstract Aligner class cannot spawn jobs.")

	def _add_postprocessor(self, result_file, data_type, command, workdir_path, job_id, delete_temps, desc_str, stdout=None):
		del_dirs= []
		if delete_temps:
			del_dirs = [workdir_path]

		postp = lambda: read_internal_alignment(result_file, datatype=data_type, del_dirs=del_dirs, file_manager=self.file_manager)

		return Job(command,
			   post_processor=postp,
                           cwd=workdir_path,
                           id=job_id,
                           description="%s_%s"%(desc_str, self.name),
                           stdout=stdout) 

class Pagan(Aligner):
	name = "pagan"
	result_suffix = ".fas"

	def __init__(self, file_manager, **kwargs):
		Aligner.__init__(self, file_manager, **kwargs)
	
	def create_job(self, alignment, guide_tree=None, **kwargs):
		unaligned_seqs = alignment.unaligned()
		job_id = kwargs.get("id", "")

		is_fake, job = self.check_and_create_fake_job(unaligned_seqs, **kwargs)
		if is_fake:
			return job

		wdir, ifn, ofn = self._preprocess_input(alignment, **kwargs)
		command = [self.cmd, "--seqfile", ifn, "--outfile", ofn]
		if guide_tree is not None:
			tree_file = os.path.join(wdir, "guide_tree.tre")
			write_tree_to_file(guide_tree, tree_file)
			command.extend(["--treefile", tree_file])
		
		if self.user_config is not None:
			command.extend(self.user_config)	

		_LOG.debug("command:%s"%" ".join(command))
		
		return self._add_postprocessor(ofn+self.result_suffix,
                                               unaligned_seqs.datatype,
                                               command,
                                               wdir,
                                               job_id,
                                               kwargs.get("delete_temps", self.delete_temps),
                                               kwargs.get("description", "whole"))


class Prank(Aligner):
	is_bundled = True
	name = "prank"
	def __init__(self, file_manager, **kwargs):
		Aligner.__init__(self, file_manager, **kwargs)
		

	def create_job(self, alignment, guide_tree=None, **kwargs):
		result_suffix = ".best.fas"
		unaligned_seqs = alignment.unaligned()
		num_taxa = unaligned_seqs.get_num_taxa()

		is_fake, job = self.check_and_create_fake_job(unaligned_seqs, **kwargs)
                if is_fake:
                        return job

		wdir, ifn, ofn = self._preprocess_input(alignment, **kwargs)	
		command = [self.cmd, "-d=%s"%ifn, "-o=%s"%ofn]  
                if "iterate" in kwargs:
                        command.append("-iterate=%d"%kwargs.get("iterate"))
                else:
                        command.append("-once")

		if guide_tree is not None:
			tree_file = os.path.join(wdir, "guide_tree.tre")
			write_tree_to_file(guide_tree, tree_file)
			command.append("-t=%s"%tree_file)
			old_tree = kwargs.get("old_tree", None)
			if old_tree is not None:
				old_tree_file = os.path.join(wdir, "old_tree.tre")
				write_tree_to_file(old_tree, old_tree_file)
				command.append("-ot=%s"%old_tree_file)
				result_suffix = ".fas"

		if self.user_config is not None:
			command.extend(self.user_config)

		if "output_options" in kwargs:
			command.extend(kwargs.get("output_options"))	

		if self.user_config is None or not "-seed" in "".join(self.user_config):
			if DEBUG:
                        	command.append("-seed=1111111")

		if self.user_config is None or not "-dnafreqs" in "".join(self.user_config):
			if alignment.datatype == "dna":
				command.append("-dnafreqs=%s"%",".join([str(freq) for freq in alignment.dna_freqs]))

		command.append("-tmp=%s"%wdir)
		if num_taxa > 200:
			command.append("-uselogs")

		_LOG.debug("PRANK:command wrapped %s"%" ".join(command))
		_LOG.debug("Prank delete_temps:%s"%kwargs.get("delete_temps", self.delete_temps))

		return self._add_postprocessor(ofn+result_suffix,
					       unaligned_seqs.datatype,
                                               command,
                                               wdir,
                                               kwargs.get("id", ""),
                                               kwargs.get("delete_temps", self.delete_temps),
                                               kwargs.get("description", "whole"))	 

	def translate(self, input_path, wdir, args=""):
		if not os.path.exists(wdir):
			os.makedir(wdir)

		file_suffix = "_translated"
		file_prefix = os.path.splitext(os.path.basename(input_path))[0] + file_suffix

		if self.user_config is not None:
			args = "%s %s"%(" ".join(self.user_config), args)

		commands="%s -d=%s -o=%s/%s -convert -translate -keep %s>>%s/stdout.txt"%(self.cmd, input_path, wdir, file_prefix, args, wdir)
		stdout_file_name = "%s/stdout.txt"%wdir
		res = os.system(commands)

		err_msg = "".join(filter(lambda x: "Unknown" in x or "Warning" in x, open(stdout_file_name, 'r')))
		if not "Warning" in err_msg and res:
                        raise RuntimeError("Error happened during PRANK translation. ErrMsg:%s"%err_msg)
                else:   
                        os.remove(stdout_file_name)

                if err_msg:
                        _LOG.info("During translating. Msg:%s, Input File:%s"%(err_msg, input_path))

		return os.path.join(wdir, filter(lambda x: file_prefix in x, os.listdir(wdir))[0])
			
		
class Mafft(Aligner):
	is_bundled = True
	name = "mafft"
	accurate_mode = True

	def __init__(self, file_manager, **kwargs):
		Aligner.__init__(self, file_manager, **kwargs)
		self.accurate_mode = (not kwargs.get('nolocal', False))

	def create_job(self, alignment, guide_tree=None, **kwargs):
		unaligned_seqs = alignment.unaligned()

		is_fake, job = self.check_and_create_fake_job(unaligned_seqs, **kwargs)
                if is_fake:
                        return job

		wdir, ifn, ofn = self._preprocess_input(alignment, **kwargs)
		command = [self.cmd]

		if self.accurate_mode and unaligned_seqs.get_num_taxa() < 200 and unaligned_seqs.max_sequence_length() < 2000:
			command.extend(["--localpair", "--maxiterate", "1000"])

		command.extend(["--quiet", ifn])
		if self.user_config is not None:
			command.extend(self.user_config)

		return self._add_postprocessor(ofn,
                                               unaligned_seqs.datatype,
                                               command,
                                               wdir,
                                               kwargs.get("id", ""),
                                               kwargs.get("delete_temps", self.delete_temps),
                                               kwargs.get("description", "whole"),
                                               ofn)

class ClustalW(Aligner):
        is_bundled = True
	name = "clustalw"

        def __init__(self, file_manager, **kwargs):
                Aligner.__init__(self, file_manager, **kwargs)

        def create_job(self, alignment, guide_tree=None, **kwargs):
                unaligned_seqs = alignment.unaligned()

		is_fake, job = self.check_and_create_fake_job(unaligned_seqs, **kwargs)
                if is_fake:
                        return job

                wdir, ifn, ofn = self._preprocess_input(alignment, **kwargs)
                command = [self.cmd, "-align", "-infile=%s"%ifn, "-outfile=%s"%ofn, "-output=fasta"]

                if self.user_config is not None:
			command.extend(self.user_config)

                return self._add_postprocessor(ofn,
                                               unaligned_seqs.datatype,
                                               command,
                                               wdir,
                                               kwargs.get("id", ""),
                                               kwargs.get("delete_temps", self.delete_temps),
                                               kwargs.get("description", "whole"))

class Merger(Tool):
	result_prefix = "merged"

	def __init__(self, file_manager, **kwargs):
		self.name += "merger"
		Tool.__init__(self, file_manager, **kwargs)

	def _preprocess_input(self, alignment1, alignment2, guide_tree1, guide_tree2, **kwargs):
		workdir_path = self.make_workdir(kwargs.get("tmp_dir"), kwargs.get("id",""))
		input_file1 = os.path.join(workdir_path, "1.fasta")
		input_file2 = None
		tree_file1 = None
		tree_file2 = None
		tree_file1 = os.path.join(workdir_path, "guide_tree1.tre")
		is_merge = kwargs.get("is_merge", True)		

		if not is_merge:
                	input_file_obj = open(input_file1, 'w')  
                	alignment1.write_to_stream(input_file_obj, "fasta", " group_A")
                	alignment2.write_to_stream(input_file_obj, "fasta", " group_B")
                	input_file_obj.close()
			write_tree_to_file(kwargs.get("tree"), tree_file1)	
		else:
			input_file2 = os.path.join(workdir_path, "2.fasta")
                        alignment1.write_to_path(input_file1)
                        alignment2.write_to_path(input_file2)	
			if guide_tree1 is not None and guide_tree2 is not None:
				tree_file2 = os.path.join(workdir_path, "guide_tree2.tre")
				write_tree_to_file(guide_tree1, tree_file1)
				write_tree_to_file(guide_tree2, tree_file2)

		output_file = os.path.join(workdir_path, self.result_prefix)

		return workdir_path, input_file1, input_file2, tree_file1, tree_file2, output_file	

	def check_and_create_fake_job(self, alignment1, alignment2, **kwargs):
		num_taxa1 = alignment1.get_num_taxa()
                num_taxa2 = alignment2.get_num_taxa()

                if num_taxa1 < 1 and num_taxa2 < 1:
                        raise RuntimeError("Neither alignment1 and alignment2 shouldn't be empty.")

                if not(num_taxa1 > 0 and num_taxa2 > 0):
                        kwargs["description"] = "%s_fake%s"%(kwargs.get("description", "whole"), self.name)
                        if num_taxa1 < 1:
                                return True, FakeJob(alignment2, **kwargs)
                        else:
                                return True, FakeJob(alignment1, **kwargs)

		return False, None
	


	def _add_postprocessor(self, result_file, datatype, command, workdir_path, job_id, delete_temps, desc_str):
		del_dirs = []
		if delete_temps:
			del_dirs.append(workdir_path)

		postp = lambda: read_internal_alignment(result_file, datatype=datatype, del_dirs=del_dirs, file_manager=self.file_manager)

		return Job(command,
			   post_processor=postp,
                           cwd=workdir_path,
                           id=job_id,
                           description="%s_%s"%(desc_str, self.name))

class MuscleMerger(Merger):
	is_bundled = True
	name = "muscle"

	def __init__(self, file_manager, **kwargs):
		Merger.__init__(self, file_manager, **kwargs)

	def create_job(self, alignment1, alignment2, guide_tree1=None, guide_tree2=None, **kwargs):
		is_fake, job = self.check_and_create_fake_job(alignment1, alignment2, **kwargs)
		if is_fake:
			return job

		wdir, input1, input2, tree1, tree2, output = self._preprocess_input(alignment1, alignment2, guide_tree1, guide_tree2, **kwargs)
		command = [self.cmd, "-in1", input1, "-in2", input2, "-out", output, "-quiet", "-profile"]
		return self._add_postprocessor(result_file=output,
					      datatype=alignment1.datatype,
                                              command=command,
                                              workdir_path=wdir,
                                              job_id=kwargs.get("id", ""),
                                              delete_temps=kwargs.get("delete_temps", self.delete_temps),
                                              desc_str=kwargs.get("description", "whole"))

class PrankMerger(Merger):
	is_bundled = True
	name = "prank"
	result_suffix = ".fas"

	def __init__(self, file_manager, **kwargs):
		Merger.__init__(self, file_manager, **kwargs)
	
	def create_job(self, alignment1, alignment2, guide_tree1=None, guide_tree2=None, **kwargs):
		is_merge = True #True: using merge function of PRANK, False:using partalin function of PRANK
                num_taxa1 = alignment1.get_num_taxa()
                num_taxa2 = alignment2.get_num_taxa()

		is_fake, job = self.check_and_create_fake_job(alignment1, alignment2, **kwargs)
		if is_fake:
			return job

		partalign_max_num = 3
		if num_taxa1 < partalign_max_num or num_taxa2 < partalign_max_num:
			is_merge = False	

                wdir, input1, input2, tree1, tree2, output = self._preprocess_input(alignment1, alignment2, guide_tree1, guide_tree2, is_merge=is_merge, **kwargs)
	
		command = [self.cmd, "-o=%s"%output]
		if is_merge:
			command.extend(["-d1=%s"%input1,"-d2=%s"%input2])
			if tree1 is not None and tree2 is not None:
				command.extend(["-t1=%s"%tree1, "-t2=%s"%tree2])
		else:
			command.extend(["-d=%s"%input1, "-partaligned"])
			if tree1 is not None:
				command.append("-t=%s"%tree1)

		if num_taxa1 + num_taxa2 > 200:
			command.append("-uselogs")

		if is_merge and kwargs.get("merge_tree", None):
			merge_tree = os.path.join(wdir, "merge_tree")
			write_tree_to_file(kwargs.get("merge_tree"), merge_tree)
			command.append("-t=%s"%merge_tree)

		if self.user_config is not None:
			command.extend(self.user_config)

		if self.user_config is None or not "-seed" in "".join(self.user_config):
                        if DEBUG:
                                command.append("-seed=1111111")

		if self.user_config is None or not "-dnafreqs" in "".join(self.user_config):
			if alignment1.datatype == "dna":
				command.append("-dnafreqs=%s"%",".join([str(freq) for freq in alignment1.dna_freqs]))

		if "output_options" in kwargs:
			command.extend(kwargs.get("output_options"))	

		command.append("-tmp=%s"%wdir)
		_LOG.debug("command:%s"%" ".join(command))	
		#Prank version 130129 rename the output name without index number
		return self._add_postprocessor(result_file=output+self.result_suffix,
				               datatype=alignment1.datatype,
                                               command=command,
                                               workdir_path=wdir,
                                               job_id=kwargs.get("id", ""),
                                               delete_temps=kwargs.get("delete_temps", self.delete_temps), 
                                               desc_str=kwargs.get("description", "whole"))
		
		
	
class TreeEstimator(Tool):
	def __init__(self, file_manager, **kwargs):
		Tool.__init__(self, file_manager, **kwargs)
				
	def _preprocess_input(self, alignment, **kwargs):
		raise NotImplementedError("Abstract TreeEstimator class.")

	@staticmethod
	def _read_results(result_file):
		raise NotImplementedError("Abstract TreeEstimator class.")

class PrankTree(TreeEstimator):
        name = "pranktree"

        def __init__(self, file_manager, **kwargs):
                TreeEstimator.__init__(self, file_manager, **kwargs)

        def _preprocess_input(self, alignment, **kwargs):
                workdir_path = self.make_workdir(kwargs.get("tmp_dir"), kwargs.get("id", ""))
                input_file = os.path.join(workdir_path, "input.fas")
                alignment.write_to_path(input_file)

                return workdir_path, input_file

        def create_job(self, alignment, **kwargs):
                wdir, input_file = self._preprocess_input(alignment, **kwargs)
                result_file = os.path.join(wdir,"njtree")
                command = [self.cmd, "-d=%s"%input_file, "-o=%s"%result_file, "-njtree", "-treeonly"]
                _LOG.debug("command:%s"%" ".join(command))

                if self.user_config is not None:
                        command.extend(self.user_config)

                def read_result(result_file):
                        return None, file(result_file).read()

                postp = lambda: read_result(result_file+".dnd")
                return Job(command,
                           post_processor=postp,
                           cwd=wdir,
                           id=kwargs.get("id",""),
                           description="%s_%s"%(kwargs.get("description", "whole"), self.name))

class FastTree(TreeEstimator):
	name = "fasttree"

	def __init__(self, file_manager, **kwargs):
		TreeEstimator.__init__(self, file_manager, **kwargs)

	def _preprocess_input(self, alignment, **kwargs):
		workdir_path = self.make_workdir(kwargs.get('tmp_dir'), kwargs.get("id", ""))
		input_file = os.path.join(workdir_path, "input.phy")
		alignment.write_to_path(input_file, "phylip")
		model = None	

		if "dna" == alignment.datatype:
			model = "-nt -gtr -gamma"
		elif "protein" == alignment.datatype:
			model = "-wag"
		else:
			raise ValueError("Datatype '%s' not suppported by FastTree"%str(alignment.datatype))	

		return workdir_path, input_file, model
	
	def create_job(self, alignment, guide_tree=None, **kwargs):
		wdir, input_file, model = self._preprocess_input(alignment, **kwargs)
		command = [self.cmd]
		command.extend(model.split())
		command.append(input_file)

		if self.user_config is not None:
			command.extend(self.user_config)


		if guide_tree is not None:
			tree_file = os.path.join(os.path.abspath(wdir), "start.tre")
			write_tree_to_file(guide_tree, tree_file)	
			command.extend(["-intree", tree_file])

		num_cpus = kwargs.get("num_cpus", 1)

		if num_cpus > 1:
			old_cmd = self.cmd
			if not self.cmd.endswith("MP"):
				self.cmd += "MP"

			res, msg = self.check_executable()
			if res:
				command[0] = self.cmd
				os.putenv("OMP_NUM_THREADS", str(num_cpus))
			else:
				self.cmd = old_cmd

		result_file = wdir+"/result.tree"

		del_dirs = []
		if kwargs.get("delete_temps", self.delete_temps):
			del_dirs = [wdir]		

		postp = lambda: read_fasttree_results(result_file, wdir, del_dirs, self.file_manager)
		_LOG.debug("command:%s"%" ".join(command))
		return Job(command,
			   post_processor=postp,
                           cwd=wdir,
                           id=kwargs.get("id", ""),
                           description="%s_%s"%(kwargs.get("description", "whole"), self.name),
                           stdout=result_file)


class Raxml(TreeEstimator):
	name = "raxml"

	def __init__(self, file_manager, **kwargs):
		TreeEstimator.__init__(self, file_manager, **kwargs)

	def _preprocess_input(self, alignment, **kwargs):
		self.name = self.name + kwargs.get("name_suffix", "")
		workdir_path = self.make_workdir(kwargs.get('tmp_dir'), kwargs.get("id", ""))
		input_file = os.path.join(workdir_path, "input.phy")
		alignment.write_to_path(input_file, "phylip")

		model = kwargs.get("model", None)
		if "dna" == alignment.datatype:
			if model is None or "" == model:
				model = "GTRGAMMA"
		elif "protein" == alignment.datatype:
			if model is None or "" == model:
				model = "PROTGAMMAWAG"
		else:
			raise ValueError("Datatype '%s' not suppported by RAxML"%str(alignment.datatype))	

		return workdir_path, input_file, model


	def create_job(self, alignment, guide_tree=None, result_suffix="default", **kwargs):
		wdir, input_file, model = self._preprocess_input(alignment, **kwargs)
		command = [self.cmd, "-m", model, "-n", result_suffix, "-s", input_file]

		if self.user_config is not None:
			command.extend(self.user_config)

		if self.user_config is None or "-p" not in self.user_config:
                        if DEBUG:
                                command.extend(["-p", "11111111"])
                        else:
                                command.extend(["-p", str(random.randint(1, int(time.time())))])


		if guide_tree is not None:
			tree_file = os.path.join(os.path.abspath(wdir), "start.tre")
			write_tree_to_file(guide_tree, tree_file)	
			command.extend(["-t", tree_file])

		num_cpus = kwargs.get("num_cpus", 1)
		if num_cpus > 1:
			old_cmd = self.cmd
			command.extend(["-T", str(num_cpus)])
			if not self.cmd.endswith("p"):
				self.cmd += "p"

			res, msg = self.check_executable()
			if res:
				command[0] = self.cmd
			self.cmd = old_cmd

		
		del_dirs = []
		if kwargs.get("delete_temps", self.delete_temps):
			del_dirs = [wdir]		

		postp = lambda: read_raxml_results(wdir, del_dirs, self.file_manager)
		_LOG.debug("command:%s"%" ".join(command))
		return Job(command,
			   post_processor=postp,
                           cwd=wdir,
                           id=kwargs.get("id", ""),
                           description="%s_%s"%(kwargs.get("description","whole"),
                           self.name))
		
				
class PhyML(TreeEstimator):
	name = "phyml"

        def __init__(self, file_manager, **kwargs):
                TreeEstimator.__init__(self, file_manager, **kwargs)

        def _preprocess_input(self, alignment, **kwargs):
                workdir_path = self.make_workdir(kwargs.get("tmp_dir"), kwargs.get("id", ""))
                input_file = os.path.join(workdir_path, "input.phy")
                data_type = "nt"
                alignment.write_to_path(input_file, "phylip")

                model = kwargs.get("model", None)
                if "dna" == alignment.datatype:
			if model is None or"" ==  model:
                        	model = "HKY85"
                        elif model not in ["HKY85", "JC69", "K80", "F81", "F84", "TN93", "GTR", "custom"]:
                                raise ValueError("The %s model is not available for PhyML."%model)
                elif "protein" == alignment.datatype:
                        data_type = "aa"
			if model is None or"" ==  model:
                        	model = "WAG"
                        elif model not in ["LG", "WAG", "JTT", "MtREV", "Dayhoff", "DCMut", "RtREV", "CpREV", "VT", "Blosum63", "MtMam", "MtArt", "HIVw", "HIVb", "custom"]:
                                raise ValueError("The %s model is not available for PhyML."%model)
                else:  
                        raise ValueError("Datatype '%s' not suppported by %s"%(str(alignment.datatype), self.name))

                return workdir_path, input_file, model, data_type

        def create_job(self, alignment, guide_tree=None, result_suffix="default", **kwargs):
                wdir, input_file, model, data_type = self._preprocess_input(alignment, **kwargs)
                command = [self.cmd, "-i", input_file, "-m", model, "-d", data_type, "--run_id", result_suffix]
                if self.user_config is not None:
                        command.extend(self.user_config)

                if self.user_config is not None and "--r_seed" not in self.user_config:
                        if DEBUG:
                                command.extend(["--r_seed", "11111111"])
                        else:  
                                command.extend(["--r_seed", str(random.randint(1, int(time.time())))])


                if guide_tree is not None:
                        tree_file = os.path.join(os.path.abspath(wdir), "start.tre")
                        write_tree_to_file(guide_tree, tree_file)
                        command.extend(["-u", tree_file])

                del_dirs = []
                if kwargs.get("delete_temps", self.delete_temps):
                        del_dirs = [wdir]

                postp = lambda: read_phyml_results(wdir, del_dirs, self.file_manager)
                _LOG.debug("command:%s"%" ".join(command))

                return Job(command,
			   post_processor=postp,
                           cwd=wdir,
                           id=kwargs.get("id", ""),
                           description="%s_%s"%(kwargs.get("description", "Whole"), self.name))

ALIGNER_CLASS = [Prank, Mafft, ClustalW, Pagan]
MERGER_CLASS = [MuscleMerger, PrankMerger]
TREE_ESTIMATOR_CLASS = [FastTree, Raxml, PhyML, PrankTree]
TRANSLATOR_CLASS= [Prank]
