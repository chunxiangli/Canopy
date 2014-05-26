import os, platform

DEBUG = False
# global debugging flag
if "DEBUG" in os.environ:
    if os.environ["DEBUG"]:
        if os.environ["DEBUG"].lower()[0] in ["1", "t", "y", "d"]:
            DEBUG = True
        else:
            DEBUG = False
    else:
        DEBUG = False
else:
    DEBUG = False

def get_number_of_cpus():
	try:
		import multiprocessing
		return multiprocessing.cpu_count()
	except (ImportError, NotImplementedError):
		pass

	#windows
	try:
		return int(os.environ("NUMBER_OF_PROCESSORS"))
	except (KeyError, ValueError):
		pass
	
	#linux
	try:
		with open("/proc/cpuinfo", 'r') as f:
			return f.read().count("processor")
	except IOError:
		pass

	raise Exception("Can't read the number of processors!")

def deploy_tools():
	if platform.system() == "Darwin":
		return os.path.realpath("iprank/mac_bin/")
	elif platform.system() == "Linux":
        	return os.path.realpath("iprank/bin/")

class Config(object):
	sub_num = 0
	pass

class Option(object):
        name = None
        value = None
        def __init__(self, short_name, long_name, type=None, group=None, **kwargs):
                assert short_name, "At least specify a short option string."

                self._short_str = short_name
		self._long_str = long_name
                self._type = type
                self._group = group
                self._kwargs = kwargs

        def add_to_parser(self, parser):
                if self._type is not None:
                        parser.add_option(self._long_str, type=self._type, **self._kwargs)
                else:
                        parser.add_option(self._long_str, **self._kwargs)

        def get_value_from_config(self, config):
                func_name = config.get
                if self._type == 'int':
                        func_name = config.getint

                if self._type == 'float':
                        func_name = config.getfloat

                if "action" in self._kwargs:
			if self._kwargs["action"] in ["store_true", "store_false"]:
                        	func_name = config.getboolean

                try:
			name = self._short_str
			if "." in name:
				name = name.split(".")[1]
                        value = func_name(self._group, name)
			if 'choice' == self._type:
                                if not value in self._kwargs.get('choices'):
                                        raise ValueError('%s value is invalid, if not sure, please check help.'%self._short_str)
                except ValueError, e:
                        raise e
                except:
                        return False
                else:
			if value != '':
                        	return {self._short_str:value}
			else:
				return False

			
class ConfigAndOptionParser(dict, object):
        _sections = dict()
        _config = None

        def __init__(self):  
                dict.__init__(self)
                self.set_sections()
	#TODO:allow multiple config files
        def read_from_config(self, config_files):
                import ConfigParser
                self._config = ConfigParser.ConfigParser()
                self._config.read(config_files)
                self.set_default_path(self._sections.keys())
                for section in self._sections.keys():
                        if self._config.has_section(section):
                                if section not in self:
                                        self[section] = {}
                                for op in self._sections[section]:
                                        d = op.get_value_from_config(self._config)
                                        if d is not False:
                                                self[section].update(d)

        def read_from_commandline(self):
                import optparse
                usage='Usage: %prog [options] file1 file2'
                option = optparse.OptionParser(usage=usage)
                self.set_options_to_parser(option)
                self._command, cfile_list = option.parse_args()
                self.read_from_config(cfile_list)
                self.update_from_dict(option, self._command.__dict__)

        def set_sections(self):
                g = []
                g.append(Option('seq', '--seq', group='main', help='Sequences file path. If the sequences file is in NEXUS format, the start tree can also be given in it.'))
                g.append(Option('work_directory', '--work_directory', group='main', help='Directory for output files and temporary files'))
 		g.append(Option('archive_directory', '--archive_directory', group='main', help='Remote directory for archived result.'))
                g.append(Option('name', '--name', group='main', help='Job name which will be used as the workdirectory name'))
                g.append(Option('datatype', '--datatype', group='main', help='The type of the input sequences'))
                g.append(Option('format', '--format', group='main', help='Input sequence files format.'))
                g.append(Option('output_format', '--output_format', group='main', type='choice', choices=['phylip', 'fasta', 'nexus'], help="Output sequence files format"))
		g.append(Option('tree_format', '--tree_format', group='main', help='Result tree file format'))
                g.append(Option('align_datatype', '--align_datatype', group='main', type="choice", choices=['PROTEIN', 'CODON', 'DNA'], help="The expected type of alignment:PROTEIN, CODON or DNA.Don't need to specify if only use the origin input sequence"))
                g.append(Option('max_iter', '--max_iter', type='int', group='main', help='The maximum number of iterations'))
                g.append(Option('max_time', '--max_time', type='int', group='main', help='The maximum time span for iterations'))
                g.append(Option('max_prob_size', '--max_prob_size', type='int', group='main', help='The maximum size of subalignment'))
                g.append(Option('num_cpus', '--num_cpus', type='int', group='main', help='The maximum number of cpus'))
                g.append(Option('need_sub_iter', '--need_sub_iter', group='main', help='Whether need iteration for subalignment', action='store_true'))
                g.append(Option('test', '--test', group='main', help='Whether need to measure the alignment and infered tree with the prior to true alignment and tree', action='store_true'))
                g.append(Option('tree', '--tree', group='main', help='Start tree file path.'))
                g.append(Option('tree_str', '--tree_str', group='main', help='Start tree string in newick format'))
                g.append(Option('score', '--score', type='float', group='main', help='Likelihood score of start tree'))
                g.append(Option('initial_aligner', '--initial_aligner', group='main', help='The name of initial aligner tool. Currently support: mafft, prank, clustalw. Dfault:mafft'))
                g.append(Option('aligner', '--aligner', group='main', help='The name of aligner tool for iteration steps. Currently support: mafft, prank, clustalw. Dfault:mafft'))
                g.append(Option('merger', '--merger', group='main', help='The name of merege tool for iteration steps. Currently support: muscle, prank. Dfault:prank'))
                g.append(Option('tree_estimator', '--tree_estimator', group='main', help='The name of tree estimate tool. Currently support: raxml, fasttree, phyml. Dfault:raxml'))
		g.append(Option('config_export_path', '--config_export_path', group='main', help="the file for saving user's configure"))
		g.append(Option('without_guide', '--without_guide', group='main', help="Don't provide guide tree to align procedure.", action='store_true'))
		g.append(Option('decomposite_strategy', '--decomposite_strategy', type='choice', choices=['longest', 'seed'], group='main', help='The strategy for the decomposition of root tree. longest: remove the longest branch to split a tree into two parts; seed: split a tree by removing the seed node.'))
		g.append(Option('save_option', '--save_option', group='main', type='choice', choices=['simple', 'rich', 'all'], help='Options for saving temp files. simple: only store input, initial and final alignments and trees; rich: additionally store alignments and trees for all iterations; all: all tmp files will be compressed in a single file.Default:simple'))
		g.append(Option('showanc', '--showanc', group='main', help='output ancestral sequences.', action='store_true'))
		g.append(Option('showevents', '--showevents', group='main', help='output evolutioanry events', action='store_true'))
		g.append(Option('showxml', '--showxml', group='main', help='output xml-files', action='store_true'))
		g.append(Option('showall', '--showall', group='main', help='output all files', action='store_true'))

                self._sections['main'] = g
                g = []
                g.append(Option('mafft.path', '--mafft.path', group='mafft', help='Prank exectutable path'))
                g.append(Option('mafft.args', '--mafft.args', group='mafft', help='Extra arguments for mafft'))
                self._sections['mafft'] = g
		g = []
                g.append(Option('clustalw.path', '--clustalw.path', group='clustalw', help='Prank exectutable path'))
                g.append(Option('clustalw.args', '--clustalw.args', group='clustalw', help='Extra arguments for clustalw'))
                self._sections['clustalw'] = g

                g = []
                g.append(Option('raxml.model', '--model', group='raxml', help='Substitution model. The default model for DNA sequences is GTRGAMMA and PROTGAMMAWAG for amino acid sequences.'))
                g.append(Option('raxml.path', '--raxml.path', group='raxml', help='Raxml exectutable path'))
                g.append(Option('raxml.args', '--raxml.args', group='raxml', help='Extra arguments'))
                self._sections['raxml'] = g

		g = []
                g.append(Option('phyml.model', '--phyml.model', group='phyml', help='Substitution model. The default model for DNA sequences is HKY85 and WAG for amino acid sequences.'))
                g.append(Option('phyml.path', '--phyml.path', group='phyml', help='phyml exectutable path'))
                g.append(Option('phyml.args', '--phyml.args', group='phyml', help='Extra arguments'))
                self._sections['phyml'] = g
		
                g = []
                g.append(Option('prank.path', '--prank.path', group='prank', help='Prank exectutable path'))
                g.append(Option('prank.args', '--prank.args', group='prank', help='Extra arguments'))
                self._sections['prank'] = g
		g = []
                g.append(Option('pagan.path', '--pagan.path', group='pagan', help='Pagan exectutable path'))
                g.append(Option('pagan.args', '--pagan.args', group='pagan', help='Extra arguments'))
                self._sections['pagan'] = g

		g = []
		g.append(Option('muscle.path', '--muscle.path', group='muscle', help='muscle exectutable path'))
                g.append(Option('muscle.args', '--muscle.args', group='muscle', help='Extra arguments'))
                self._sections['muscle'] = g

		g = []
                g.append(Option('pranktree.path', '--pranktree.path', group='pranktree', help='pranktree exectutable path'))
                g.append(Option('pranktree.args', '--pranktree.args', group='pranktree', help='Extra arguments'))
                self._sections['pranktree'] = g

		g = []
		g.append(Option('fasttree.path', '--fasttree.path', group='fasttree', help='fasttree exectutable path'))
                g.append(Option('fasttree.args', '--fasttree.args', group='fasttree', help='Extra arguments'))
                self._sections['fasttree'] = g




        def set_options_to_parser(self, parser):
                from optparse import OptionGroup
                for section in self._sections.keys():
                        group = OptionGroup(parser, section, 'options for '+section)
                        for p in self._sections[section]:
                                p.add_to_parser(group)
                        parser.add_option_group(group)
	
	def update_from_dict(self, parser, dict_option):
                for k in dict_option.keys():
                        v = dict_option[k]
                        if v is not None:
                                long_str = '--'+k
                                section = parser.get_option_group(long_str).title

                                if not section in self:
                                        self[section] = {}
                                self[section].update({k: v})

        def set_default_path(self, sections):
                bin_dir = deploy_tools()
                for s in sections:
                        if s != 'main':
                                self[s] = {'%s.path'%s: os.path.join(bin_dir, s)}
                                if platform.system() == 'Windows':
                                        self[s].path += '.exe'

        def write_to_file(self, file_path):
		file_path_dir = os.path.dirname(file_path)
		if "" != file_path_dir and os.path.exists(os.path.realpath(file_path_dir)):
                	file_real_path = os.path.realpath(file_path)
		else:
			file_real_path = os.path.join(self["main"]["work_directory"], os.path.basename(file_path))

		self["main"]["config_export_path"] = file_real_path
                file_obj = open(file_real_path, 'w')

		for k, v in self.iteritems():
			if isinstance(v,dict):
				file_obj.write("[%s]\n"%k)
				for k1, v1 in v.iteritems():
					if v1 is not None:
						k1_arr = k1.split('.')
						if len(k1_arr) > 1:
							k1= k1_arr[1]
						file_obj.write("%s=%s\n"%(k1,v1))
				file_obj.write("\n")
                file_obj.close()

		return file_real_path
