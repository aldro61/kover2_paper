import argparse
import os
import subprocess
import numpy as np
import warnings

import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 8}
matplotlib.rc('font', **font)

class TreeNode(object):
	def __init__(self, kmer=None):
		self.kmer = kmer
		self.parent = None
		
	@property
	def is_leaf(self):
		return self.kmer is None
		
class ModelParser(object):
	
	def __init__(self, input_path):
		self.input_path = input_path
		self.tree = None
		self.parse()
		
	def get_tree(self):
		if self.tree is None:
			raise RuntimeError("Tree is None")
		else:
			return self.tree
		
	def extract_rule_id(self, header):
		return int(header.split(",")[0].split(":")[-1].split("___")[0])
		
	def extract_rule_node(self, header, kmer):
		spt = header.split(",")
		ex = spt[0].split(":")[-1].split("___")[1].split("_")[1]
		eq = spt[0].split(":")[-1].split("___")[2].split("_")[1]
		imp = float(spt[-1].split(":")[-1].lstrip())
		node = TreeNode(kmer=kmer)
		node.nb_examples = ex
		node.eq_rules = eq
		node.importance = imp
		return node
		
	def extract_class_probas(self, header):
		spt = header.split("___")[2].split("__")
		classes, probas = zip(*[(class_info.split("_")[0], float(class_info.split("_")[2])) for class_info in spt])
		probas = np.array(probas)
		rounded_probas = probas.round(2)
		error = 1.0 - np.sum(rounded_probas)
		n = int(round(error / 0.01))
		fixed_probas = rounded_probas[:]
		for _,j in sorted(((probas[j] - rounded_probas[j], j) for j in range(probas.shape[0])), reverse=n>0)[:abs(n)]:
			fixed_probas[j] += math.copysign(0.01, n)
		return {c:p for c,p in zip(classes, fixed_probas)}
		
	def extract_children(self, header):
		spt = header.split(",")
		left = spt[1].split(":")[-1].strip()
		right = spt[2].split(":")[-1].strip()

		try:
			left = int(left.split("___")[0])
			left = self.rule_dict[left]
		except:
			class_probas = self.extract_class_probas(left)
			ex = left.split("___")[1].split("_")[1]
			left = TreeNode()
			left.nb_examples = ex
			left.resistant = class_probas["resistant"]
			left.sensitive = class_probas["sensitive"]
		
		try:
			right = int(right.split("___")[0])
			right = self.rule_dict[right]
		except:
			class_probas = self.extract_class_probas(right)
			ex = right.split("___")[1].split("_")[1]
			right = TreeNode()
			right.nb_examples = ex
			right.resistant = class_probas["resistant"]
			right.sensitive = class_probas["sensitive"]
			
		return left, right
		
	def find_root(self):
		roots = [id for id in self.rule_dict if self.rule_dict[id].parent is None]
		assert len(roots) == 1
		return roots[0]
			
	def parse(self):
		self.rules = list(zip(*fasta_to_contigs(self.input_path, return_headers=True)[::-1]))
		self.rule_dict = {self.extract_rule_id(header): self.extract_rule_node(header, kmer) for header, kmer in self.rules}
		
		for header, kmer in self.rules:
			current_rule = self.rule_dict[self.extract_rule_id(header)]
			
			left, right = self.extract_children(header)
			
			current_rule.left = left
			current_rule.right = right
			
			left.parent = current_rule
			right.parent = current_rule
		
		self.tree = self.rule_dict[self.find_root()]

def visualize_model(input_path):
	tree = ModelParser(input_path).get_tree()
	try:
		tree = ModelParser(input_path).get_tree()
	except:
		warnings.warn("Unable to parse {}. Invalid decision tree model".format(input_path))
		return
	latex_str = _latex_export(tree)
	base_filename = os.path.splitext(input_path)[0]
	tex_filename = base_filename + ".tex"
	with open(tex_filename, 'w') as tex_file:
		tex_file.write(latex_str)
	
	subprocess.call(["lualatex", os.path.basename(tex_filename)], cwd=os.path.join(".", os.path.dirname(tex_filename))) 
	os.remove(tex_filename)
	for ext in [".log", ".aux", ".synctex.gz"]:
		if os.path.exists(base_filename + ext):
			os.remove(base_filename + ext)
		
def fasta_to_contigs(path, return_headers=False):
    """
    Reads a FASTA file and loads its contigs
    Note: sequences are returned in lower case
    """
    contigs = []
    headers = []
    def add_contig(header, seq):
        if len(seq) == 0:
            raise Exception("Attempted to add a contig with length 0. Not normal! Path is " + path)
        contigs.append(seq.lower())
        headers.append(header.lower())

    buffer = ""
    header = ""
    for l in open(path, "r"):
        if l.startswith(">"):
            # New contig starting

            # Save current sequence and flush buffer
            if len(buffer) > 0.:
                add_contig(header, buffer)
                buffer = ""

            # Read contig header
            header = l[1:].strip()

        else:
            # Accumulate DNA sequence
            buffer += l.strip()

    # Save final buffer
    if buffer is not None and buffer != "":
        add_contig(header, buffer)

    if return_headers:
        return contigs, headers
    else:
        return contigs
		
def _latex_export(model):
	leaf_dict = {}
	def _rec_export(node, depth):
		if not node.is_leaf:
			indent = "\t" * depth
			return '\n{0!s} {1!s}[as=\\scriptsize {2!s}\\\\\\scriptsize {3!s}\\\\\\textbf{{{4!s}}}\\\\\\small{5!s}, nonterminal] -> {{{6!s}, {7!s}}}'\
					.format(indent,
						str(hash((node.kmer, node.parent))).replace("-", "1"),
						str(node.kmer).replace("<=", "$\leq$").replace("[", "(").replace("]", ")").replace("%", "").upper()[:16],
						str(node.kmer).replace("<=", "$\leq$").replace("[", "(").replace("]", ")").replace("%", "").upper()[16:],
						str(node.eq_rules),
						str(node.importance),
						_rec_export(node.left, depth + 1),
						_rec_export(node.right, depth + 1))
		else:
			hashed_node = str(hash(node)).replace("-", "1")
			leaf_dict[hashed_node] = ('{0!s}/green!75!black/Sensitive, {1!s}/red/Resistant' \
										.format(str(int(100*node.sensitive)), str(int(100*node.resistant))),
										str(node.nb_examples))
			return "{0!s}[as=L, terminal]".format(hashed_node)
	tree_graph = _rec_export(model, 0)
	leaf_donuts = ""
	for hashed_node, donut_params in leaf_dict.items():
		leaf_donuts += '\n\\ExtractCoordinate{{$({0!s})$}};'.format(hashed_node)
		leaf_donuts += '\n\\wheelchart{{{0!s}}}{{{1!s}}}'.format(donut_params[0], donut_params[1])
	exported = \
"""
% !TeX program = lualatex
\\RequirePackage{{luatex85}}
\\documentclass[tikz,border=5]{{standalone}}
\\usetikzlibrary{{arrows}}
\\definecolor{{nonterminal}}{{RGB}}{{230,230,230}}
\\definecolor{{terminal}}{{RGB}}{{255,51,76}}
\\usetikzlibrary{{shapes.misc, positioning}}
\\usetikzlibrary{{graphs,graphdrawing,arrows.meta}}
\\usetikzlibrary{{calc}}
\\usegdlibrary{{trees}}
\\begin{{document}}
		
\\newdimen\\XCoord
\\newdimen\\YCoord

\\newcommand*{{\\ExtractCoordinate}}[1]{{\\path (#1); \\pgfgetlastxy{{\\XCoord}}{{\\YCoord}};}}%
\\newcommand*{{\\LabelCurrentCoordinate}}[2]{{\\fill [#1] ($(\\XCoord,\\YCoord)$) circle (2pt) node [right] {{#2}}}}%
	
% Donut plot
% Adjusts the size of the wheel:
\\def\\innerradius{{0.4cm}}
\\def\\outerradius{{0.6cm}}

% The main macro
\\newcommand{{\\wheelchart}}[2]{{
	% Calculate total
	\\pgfmathsetmacro{{\\totalnum}}{{0}}
	\\foreach \\value/\\colour/\\name in {{#1}} {{
		\\pgfmathparse{{\\value+\\totalnum}}
		\\global\\let\\totalnum=\\pgfmathresult
	}}
	
	% Calculate the thickness and the middle line of the wheel
	\\pgfmathsetmacro{{\\wheelwidth}}{{\\outerradius-\\innerradius}}
	\\pgfmathsetmacro{{\\midradius}}{{(\\outerradius+\\innerradius)/2}}
	
	% Rotate so we start from the top
	\\begin{{scope}}[rotate=90]
	
	% Loop through each value set. \\cumnum keeps track of where we are in the wheel
	\\pgfmathsetmacro{{\\cumnum}}{{0}}
	\\foreach \\value/\\colour/\\name in {{#1}} {{
		\\pgfmathsetmacro{{\\newcumnum}}{{\\cumnum + \\value/\\totalnum*360}}
		
		% Calculate the percent value
		\\pgfmathsetmacro{{\\percentage}}{{\\value/\\totalnum*100}}

		% Draw the color segments. Somehow, the \\midrow units got lost, so we add 'pt' at the end. Not nice...
		\\fill[\\colour] ([shift=(-\\cumnum:\\outerradius)]\\YCoord,-\\XCoord) arc (-\\cumnum:-(\\newcumnum):\\outerradius) --
		([shift=(-\\newcumnum:\\innerradius)]\\YCoord,-\\XCoord) arc (-\\newcumnum:-(\\cumnum):\\innerradius) -- cycle;
		
		% Set the old cumulated angle to the new value
		\\global\\let\\cumnum=\\newcumnum
	}}
	
	\\end{{scope}}
	\\draw[white, fill] ($(\\XCoord,\\YCoord)$) circle (\\innerradius);
	\\draw[black] ($(\\XCoord,\\YCoord)$) node [align=center, text width=2*\\innerradius] {{#2}};
	% Uncomment for contour
	\\draw[black] ($(\\XCoord,\\YCoord)$) circle (\\outerradius) circle (\\innerradius);
}}
\\begin{{tikzpicture}}[>=Stealth,
nonterminal/.style={{draw, fill=nonterminal, rectangle, rounded corners=5pt, minimum width=1.5cm, align=center}},
terminal/.style={{draw, fill=terminal, circle, minimum size=1.2cm, align=left}}]
\\linespread{{0.75}}
\\graph[binary tree layout]{{
{0!s}
}};
{1!s}
\\end{{tikzpicture}}
\\end{{document}}
""".format(tree_graph, leaf_donuts)
	return exported

def main():
	# Parser
	parser = argparse.ArgumentParser(description="Model Interpretabiliy")
	parser.add_argument('input', type=str, default='model.fasta')
	args = parser.parse_args()
	
	visualize_model(args.input)

if __name__ == '__main__':
	main()
