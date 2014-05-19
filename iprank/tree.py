from iprank import msg_exit
dendropy_minimum_version = "3.10.0"
try:
        import dendropy
        if dendropy.__version__ < dendropy_minimum_version:
                msg_exit("Sorry, reguires Dendropy version %s or greater"%dendropy_minimum_version)
        del dendropy
except ImportError:
        msg_exit("Error, dendropy is not installed (dendropy >=%s)"%dendropy_minimum_version)

import os, copy, re
from StringIO import StringIO
from dendropy import TaxonSet, Tree, Node, treecalc, treesplit
import config

def  write_tree_to_file(tree, file_name, schema="newick"):
	with open(file_name, "w") as file_obj:
		if isinstance(tree, str):
			tree = PhylogeneticTree.read_from_string(tree)
		tree.write(file_obj, schema=schema, suppress_rooting=True)

def reroot_at_midpoint(tree_str,schema="newick"):
	phy_tree = PhylogeneticTree.read_from_string(tree_str,schema)
	phy_tree.reroot_at_midpoint(update_splits=True)
	phy_tree.resolve_polytomies()
	return phy_tree.as_newick_string()

#TEST:for sate tree
def scale_tree_branch(tree, format="newick"):
        tree_obj = None
        if os.path.exists(tree):
                tree_obj = Tree.get_from_path(tree, format)
        elif isinstance(tree, str):
                tree_obj = Tree(stream=StringIO(tree), schema=format)
        elif isinstance(tree, Tree):
                tree_obj = Tree
        if sum([ e.length > 1 for e in tree_obj.postorder_edge_iter()]):
                for e in tree_obj.postorder_edge_iter():
                        if e.length is not None:
                                e.length = e.length/100
        return tree_obj.as_newick_string()


class PhylogeneticTree(object):
	def __init__(self, dendropy_tree):
		self._tree = dendropy_tree
		self._tree.update_splits(delete_outdegree_one = False)
		self._tree.seed_node.edge.length = None
		self._tree.seed_node.edge.tail_node = None

	@property
	def tree(self):
		return self._tree

	@tree.setter
	def tree(self, t):
		self._tree = t

	@staticmethod
	def read_from_path(filename, schema="newick", taxon_set=None):
		t = Tree(taxon_set=taxon_set)
                t.read_from_path(filename, schema)

		return PhylogeneticTree(t)

	@staticmethod
	def read_from_string(tree_str, schema="newick", taxon_set=None):
		taxon = taxon_set
		if taxon is None:
			taxon = TaxonSet()

		return PhylogeneticTree(Tree.get_from_string(tree_str, schema=schema, taxon_set=taxon))

	def get_longest_internal_edge(self):
		longest_edge_length = -1
		longest_internal_edge = None

		for e in self._tree.postorder_edge_iter():
			if e.is_internal() and e.length >= longest_edge_length:
				longest_edge_length = e.length
				longest_internal_edge = e 

		return longest_internal_edge

	def bipartition_by_edge(self, edge):
		t2_len = 0
                t1_len = edge.head_node.edge.length
                t1 = None
                t2 = None

                if edge.tail_node == self._tree.seed_node:
                        edge1, edge2 =[e for e in edge.tail_node.incident_edges() if e.length != None]
                        if edge1 != edge:
                                t2_len = edge1.length
                                edge2, edge1 = edge1, edge
                        else:
                                t2_len = edge2.length
                        t2 = PhylogeneticTree(Tree(seed_node=edge2.head_node))
                        t1 = PhylogeneticTree(Tree(seed_node=edge1.head_node))
                else:
                        t = copy.deepcopy(self._tree)
                        self.prune_subtree(edge.head_node)
                        t2 = PhylogeneticTree(Tree(seed_node=edge.head_node))
                        t1 = PhylogeneticTree(self._tree)
                        self._tree = t

                return t1,t2, t1_len, t2_len


	def bipartition_by_node(self, node):
		assert not node.is_leaf(), "Cann't split a tree at a leaf node."	

		children = node.child_nodes()	
		if children[0].edge is not None:
			return self.bipartition_by_edge(children[0].edge)
		elif len(children) > 1 and children[1].edge is not None:
			return self.bipartition_by_edge(children[1].edge)
			
			
	def bipartition_by_longest_internal_edge(self):
                e = self.get_longest_internal_edge()

                assert e.head_node is not None
                assert e.tail_node is not None, self.as_newick_string()

		return self.bipartition_by_edge(e)

	def bipartition_by_seed(self):
		return self.bipartition_by_node(self._tree.seed_node)

	def prune_subtree(self, node, update_splits=False, delete_outdegree_one=True):
		self._tree.prune_subtree(node, update_splits=update_splits, delete_outdegree_one=delete_outdegree_one)

	def get_subtree_with_labels(self, labels):
                cur_leaf_labels = self.leaf_node_names()

                prune_labels = [re.sub('_',' ', label) for label in cur_leaf_labels if label not in labels]
                sub_tree = copy.deepcopy(self._tree)
                sub_tree.prune_taxa_with_labels(prune_labels)

                return PhylogeneticTree(sub_tree)


	def reroot_at_midpoint(self, update_splits=False, delete_outdegree_one=True):
		'''
			Modified from the source code of Dendropy v3.12.0.
		'''	
		pdm = treecalc.PatristicDistanceMatrix(self._tree)
		n1,n2 = pdm.max_dist_nodes
		plen = float(pdm.max_dist)/2
		mrca_node = pdm.mrca(n1.taxon, n2.taxon)
		cur_node = n1
		break_on_node = None
		target_edge = None
		head_node_edge_len = None

		while cur_node is not mrca_node:
        		if cur_node.edge.length > plen:
                		target_edge = cur_node.edge
                		head_node_edge_len = plen
                		plen = 0
                		break
        		elif abs(cur_node.edge.length - plen) < 1e-6:
                		break_on_node = cur_node.parent_node
                		break
        		else:
               			plen -= cur_node.edge.length
                		cur_node = cur_node.parent_node

		assert break_on_node is not None or target_edge is not None

		if break_on_node:
        		self._tree.reseed_at(break_on_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
		else:  
        		tail_node_edge_len = target_edge.length - head_node_edge_len
        		old_head_node = target_edge.head_node
        		old_tail_node = target_edge.tail_node
        		old_tail_node.remove_child(old_head_node)
        		new_seed_node = Node()
        		new_seed_node.add_child(old_head_node, edge_length =head_node_edge_len)
        		old_tail_node.add_child(new_seed_node, edge_length = tail_node_edge_len)
        		self._tree.reseed_at(new_seed_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)

		self._tree.is_rooted = True
	
		if update_splits:
			self._tree.update_splits(delete_outdegree_one = False)

		return self._tree.seed_node

	def leaf_node_names(self):
		return [re.sub(' ', '_', i.taxon.label) for i in self._tree.leaf_nodes()]
	
	def rename_leaf_names(self, name_map, restore=True):
                for ln in self._tree.leaf_nodes():
			if restore:
				old_name = re.sub(' ', '_', ln.taxon.label)
				if old_name in name_map:
                        		ln.taxon.label = name_map[old_name]
			else:
				ln.taxon.label = [ re.sub('_', ' ', k) for k, v in name_map.iteritems() if v == re.sub(' ', '_', ln.taxon.label)][0]

	def num_taxa(self):	
		return len(self._tree.leaf_nodes())

	def resolve_polytomies(self, update_splits=False, rng=None):
		self._tree.resolve_polytomies(update_splits, rng)

	def get_tree(self):
		return self._tree
		
	def print_plot(self):
		self._tree.print_plot()

	"debug"
	def as_newick_string(self):
		return self._tree.as_newick_string()

	def write(self, stream, schema="newick", **kwargs):
		self._tree.write(stream, schema, **kwargs)

def robinson_foulds_distance(tree_str1, tree_str2):
	if tree_str1 is None or tree_str2 is None or "" == tree_str1 or "" == tree_str2:
		return -1
	
	taxon = TaxonSet()

	return treecalc.robinson_foulds_distance(Tree(stream=StringIO(tree_str1), schema='newick', taxon_set=taxon), Tree(stream=StringIO(tree_str2), schema='newick', taxon_set=taxon))
	
	
def symmetric_difference(tree_str1, tree_str2):
	if tree_str1 is None or tree_str2 is None or "" == tree_str1 or "" == tree_str2:
		return -1

	taxon = TaxonSet()
	return treecalc.symmetric_difference(Tree(stream=StringIO(tree_str1), schema='newick', taxon_set=taxon), Tree(stream=StringIO(tree_str2), schema='newick', taxon_set=taxon))
