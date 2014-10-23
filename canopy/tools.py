from canopy import msg_exit
Bio_minimum_version = "1.58"
try:
        import Bio
        if Bio.__version__ < Bio_minimum_version:
                msg_exit("Sorry, reguires Biopython version %s or greater"%Bio_minimum_version)
        del Bio
except ImportError:
        msg_exit("Error, Biopython is not installed (Biopython >=%s)"%Bio_minimum_version)

import os, re, commands
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from canopy.logger import get_logger
from canopy.tree import PhylogeneticTree
from canopy.file_manage import remove_files

_LOG = get_logger(__name__)

INDEL_CHAR = '-'
_INDEL = re.compile(r'-')
_ILLEGAL_TREE_CHARACTER = re.compile(r'[\(\)\[\]\;]')
_DANGEROUS_CHARACTER = re.compile(r'[^#-_a-zA-Z0-9]')

def format_alignment_result_in_directory(directory, format="phylip", datatype="DNA", filter_files=None):
	file_list = os.listdir(directory)
	target_files = [ os.path.join(directory, f) for f in file_list if f.endswith(".fas") and not f.endswith("anc.fas")]
	if filter_files is not None:
		target_files = [ f for f in target_files if f not in filter_files ]
	_LOG.debug("target_files:%s"%target_files)
	file_suffix = ".phy"
	if format == "nexus":
		file_suffix = ".nex"
	name_map = {}

	for f in target_files:
		align = Alignment()
		name_map = align.read_from_path(f, data_type=datatype, name_map=name_map)

		result_file = f[:-4] + file_suffix	
		if format == "phylip":
			align.write_to_path(result_file, file_format=format) 
		else:
			align.write_to_path(result_file, file_format=format, name_map=name_map) 

		if format == "nexus":
			tree_file = f[:-4] + ".dnd"	
			if os.path.exists(tree_file):
				phy_tree = PhylogeneticTree.read_from_path(tree_file)
				with open(result_file, 'a') as rfo:
						rfo.write("begin trees;\ntree PRANK = %s\nend;"%phy_tree.as_newick_string())

	if format == "phylip" and name_map.keys() != name_map.values():
		with open(os.path.join(directory, "name_map.txt"), "w") as name_map_file:
                                name_map_file.write("\n".join(["%s %s"%(name, origin_name) for name, origin_name in name_map.iteritems()]))

	remove_files(target_files)
	
def alignment_accuracy(align_file1, align_file2):
	res = commands.getoutput("metal %s %s"%(align_file1, align_file2))
	return float(res.split("=")[1])

def alignment_accuracy_fastsp(align_file1, align_file2):
        res = commands.getoutput("java -jar ../bin/FastSP_1.3.jar -r %s -e %s"%(align_file1, align_file2))
        return float(res.split("\n")[-5].split()[1])

def name_filter_and_encode(name, name_array, name_map):
	#filter characters that crash phylogenetic tree file
	name = "".join(_ILLEGAL_TREE_CHARACTER.split(name))
	new_name = "".join(_DANGEROUS_CHARACTER.split(name))[0:10]
	nn = new_name
	for k, v in name_map.iteritems():
		if v == name:
			if k in name_array:
				raise ValueError("%s has two sequences."%name)
			else:
				name_array.append(k)
				return k

	existed_arr = [ n for n in name_map.keys() if n == nn ]
        suffix =  0
        reverse_index = 1
	while len(existed_arr):
		suffix += 1
		if suffix >= 10:
			reverse_index = len(str(suffix))
		nn = "%s%d"%(new_name[:-reverse_index], suffix)
		existed_arr = [ n for n in name_map.keys() if n == nn ]

	name_array.append(nn)
	name_map[nn] = name
	return nn

def alignment_fas(align_file):
	res  = file(align_file).readlines()
	i = 0
	name = None
	seq = ""
	a = {}
	while i < len(res):
		if res[i].startswith(">"):
			if name is not None:
				a[name] = seq
			name = res[i].split(">")[1].strip()
			seq = ""
		else:
			seq += res[i].strip()
		i += 1

def alignment_phy(align_file):
	res = file(align_file).readlines()
	num_taxa = int(res[0].split()[0])
	phylipi = False
	i = 1
	a = {}
	while i <= num_taxa:
		r = res[i].strip().split()
		name = r[0]
		a[name] = " ".join(r[1:]).replace(" ","")
		i += 1
	while i < len(res) - 1:
		phylipi = True
		i += 1
		for j in range(num_taxa):
			a[a.keys()[j]]+=res[i+j].strip().replace(' ','')
		i += num_taxa

DATATYPE_LIST = ["dna", "DNA", "protein", "PROTEIN", "rna", "RNA"]
class Alignment(dict, object):
	"""
		A simple class that maps taxa names to sequences.
	"""
	def __init__(self, align_datatype=None):
		dict.__init__(self)
		self._datatype = None
		self.names = []
		self.dna_freqs = None
		self.align_datatype = align_datatype

	@property
	def datatype(self):
		return self._datatype

	@datatype.setter	
	def datatype(self,d):
		if d is None:
			self._datatype = None
		elif d in DATATYPE_LIST:
			self._datatype = d.lower()
		else:
			raise ValueError("Datatype %s is not supported."%d)

        def reset(self):
                self.names = []
                self.clear()

	def update(self, align):
		self._datatype = align.datatype

		for name in align.names:
			if name not in self.names:
				self.names.append(name)
			self[name] = align[name]

		if self.dna_freqs is None:
			self.dna_freqs = align.dna_freqs
			
	def append(self, align):
		length = self.alignment_length()
		is_aligned = self.is_aligned()

		for name in align.names:
			try:
				self[name] += align[name]
			except KeyError:
				if is_aligned:
					self[name] = INDEL_CHAR*length+align[name] 
				else:
					self[name] = align[name]
				self.names.append(name)

		if self.dna_freqs is not None:
			self.dna_freqs = [ (self.dna_freqs[i] + align.dna_freqs[i])/2 for i in range(4)]

	def get_sequence_names(self):
		return self.names

	def get_num_taxa(self):
		return len(self.names)

	def record_name(self, name, name_map=None):
	        return	name_filter_and_encode(name, self.names, name_map)

	def read_from_path(self, filename, file_format="fasta", data_type="DNA", keys=None, name_map=None, prefix=""):
                if self.get_num_taxa() > 1:
                        self.reset()

		self.datatype= data_type

		assert os.path.exists(filename), "The sequences file %s doesn't exist."%filename 

		file_obj = open(filename, "rU")
		return self.read_from_stream(file_obj, schema=file_format, data_type=self._datatype, keys=keys, name_map=name_map)

	def read_from_stream(self, file_obj, schema="fasta", data_type="dna", keys=None, name_map=None):
		self.datatype = data_type

		if name_map is None:
			name_map = {}

		records = SeqIO.parse(file_obj, schema)

		for rec in records: 
			new_name =self.record_name(rec.id, name_map=name_map)
			if keys is None or new_name in keys:
				seq_str = str(rec.seq)

				if "CODON" == self.align_datatype and "dna" == self.datatype:
					assert not len(seq_str)%3, "The length of sequence '%s' in file '%s' is not multiple of three!"%(rec.id, os.path.basename(file_obj.name))

				self[new_name] = seq_str

		if keys is not None:
			self.names = keys

		if len(self) and self.datatype == "dna":
			self.set_dna_freqs()

		file_obj.close()
		return name_map

	def write_to_path(self, filename, file_format="fasta", name_suffix="", suppress_ancester=False, name_map=None):
		file_obj = open(filename, 'w')
		self.write_to_stream(file_obj, schema=file_format, name_suffix=name_suffix, suppress_ancester=suppress_ancester, name_map=name_map)
		file_obj.close()

	
	def write_to_stream(self, stream, schema="fasta", name_suffix="", suppress_ancester=False, name_map=None):
		records = []
		restore = False
		alphabet_dict = {"dna":generic_dna, "protein":generic_protein}

		if schema !="phylip" and isinstance(name_map, dict):
			restore = True

		def get_name(name):
			try:
				return name_map[name]
			except:
				return name

                if suppress_ancester:
                        records=[SeqRecord(Seq(self[name], alphabet_dict[self._datatype]), id=get_name(name)+name_suffix, description="") for name in self.names if not name.startswith("#")]
                else:
                        records=[SeqRecord(Seq(self[name], alphabet_dict[self._datatype]), id=get_name(name)+name_suffix, description="") for name in self.names]

		SeqIO.write(records, stream, schema)

	def unaligned(self):
		"""
			Return a new alignment with all gaps and missing sequences removed.
		"""
		new_alignment = Alignment(self.align_datatype)
		new_alignment.datatype = self._datatype
		for name in self.names:
			new_seq = re.sub(_INDEL, '', self[name])
			if new_seq != '':
				new_alignment[name] = new_seq
		new_alignment.names = self.names

		return new_alignment

	def sub_alignment(self, keys):
		sub_alignment = Alignment(self.align_datatype)
		sub_alignment.datatype = self._datatype
		sub_alignment.dna_freqs = self.dna_freqs

		for key in keys:
			if key in self.names:
				sub_alignment[key] = self[key]
				sub_alignment.names.append(key)
			else:
				_LOG.debug("Make sure the %s name is correct, which not existed in the existing alignment."%key)

		return sub_alignment

	def is_empty(self):
		return self.__len__() < 1

	def is_aligned(self):
		if self.is_empty():
			raise ValueError("The alignment is empty.\n")
		else:
			sequences = self.values()
			first_seq_len = len(sequences[0])
			return all([len(seq) == first_seq_len for seq in sequences])

	def alignment_length(self):
		if self.is_aligned():
			return len(self.values()[0])
	
	def max_sequence_length(self):
		return max(len(seq) for seq in self.values())

	def set_dna_freqs(self):
                dna_freq = [ 1 for i in xrange(4)]
                base_list = ['A', 'C', 'G', 'T']
                if self[self.names[0]][0].islower():
                        base_list = [ a.lower() for a in base_list]

                def base_count(seq):
                        def count(i):
                                dna_freq[i] = dna_freq[i] + 1

                        for i in xrange(4):
                                dna_freq[i] = dna_freq[i] + seq.count(base_list[i])
		map(base_count, self.values())

		total_num = sum(dna_freq)
		self.dna_freqs = [ (num*1.0)/total_num for num in dna_freq]

class MultiAlignments(dict, object):
	def __init__(self, align_datatype=None):
		dict.__init__(self)
		self.names = []
		self._num_taxa = 0
		self.align_datatype = align_datatype

	@property
	def datatype(self):
		return self._datatype

	@datatype.setter	
	def datatype(self,d):
		if d is None:
			self._datatype = None
		elif d in DATATYPE_LIST:
			self._datatype = d.lower()
		else:
			raise ValueError("Datatype %s is not supported."%d)

        def reset(self):
                self.names = []
                self._num_taxa = 0
                self.clear()

	def record_name(self, name):
		#return name_filter_and_encode(name, self.names, self.name_map)
		new_name = "".join(_DANGEROUS_CHARACTER.split(name))
		self.names.append(new_name)
		return new_name

	def read_from_path(self, file_list, file_format="fasta", data_type="DNA", prefix=""): 
		if self.num_taxa > 0:
                        self.reset()

		self.datatype = data_type
		name_map = {}

		for filename in file_list:
			alignment_name = os.path.splitext(os.path.basename(filename))[0].split("_translated")[0]
			if prefix != "":
				alignment_name = alignment_name.split(prefix)[1]
			new_name = self.record_name(alignment_name)
			new_alignment = Alignment(self.align_datatype)
			name_map = new_alignment.read_from_path(filename, file_format=file_format, data_type=self._datatype, name_map=name_map)
			if not new_alignment.is_empty():
				self[new_name] = new_alignment
			self._num_taxa = len(name_map)
		return name_map

	def write_to_path(self, file_path, file_format="fasta", name_suffix="", suppress_ancester=False, name_map=None):
		#TODO:Maybe only give the ouput directory
		if isinstance(file_path, list):
			num_file_path = len(file_path)
			num_alignments = len(self)

			assert num_file_path == num_alignments, "There are only %d finumames which is not enough for %d alignments."%(num_finumames, num_alignments)

			for i in range(num_file_path):
				alignment = self[self.names[i]]
				with open(file_path[i], "w") as file_stream:
					alignment.write_to_stream(file_stream,
								  schema=file_format,
							          name_suffix=name_suffix,
						                  suppress_ancester=suppress_ancester,
								  name_map=name_map)

		elif isinstance(file_path, str):
			concatenated_alignment = self.concatenate()
			with open(file_path, 'w') as file_stream:
				concatenated_alignment.write_to_stream(file_stream,
								       schema=file_format,
                                                                       name_suffix=name_suffix,
                                                                       suppress_ancester=suppress_ancester,
                                                                       name_map=name_map)

	def concatenate(self):
		concatenated_alignment = Alignment(self.align_datatype)
		concatenated_alignment.datatype = self[self.names[0]].datatype
		is_aligned = self[self.names[0]].is_aligned()

		try:
			for alignment_name in sorted(self.names):
				old_names = concatenated_alignment.names
				alignment = self[alignment_name]
				if len(concatenated_alignment) < 1:
					concatenated_alignment.update(alignment)
				else:
					concatenated_alignment.append(alignment)
				append_length = alignment.alignment_length()
				append_names = alignment.names
				
				if is_aligned:#complete the alignment with gaps
					for name in old_names:
						if not name in append_names:
							concatenated_alignment[name]+=append_length*INDEL_CHAR
		except Exception as e:
			_LOG.error("Error happend during concatenate:%s"%str(e))
			raise e
				 
		return concatenated_alignment


	def sub_alignment(self, keys):
		sub_align = MultiAlignments(self.align_datatype)
		sub_align.names = self.names

		for alignment_name in self.names:
			sub_align[alignment_name] = self[alignment_name].sub_alignment(keys)

		sub_align.num_taxa = len(keys)

		return sub_align
	
	def is_empty(self):
		for name in self.names:
			if self[name].__len__ < 1:
				return True

	@property
	def num_taxa(self):
		return self._num_taxa

	@num_taxa.setter
	def num_taxa(self, num):
		self._num_taxa = num
