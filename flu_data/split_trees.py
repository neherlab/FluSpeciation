in_tree = snakemake.input.tree
out_yam = snakemake.output.yam_all
out_vic = snakemake.output.vic

from Bio import Phylo
T = Phylo.read(in_tree, "newick")
leaves = {}
for n in T.get_terminals():
	leaves[n.name] = n
	n.tip_count=1

for n in T.get_nonterminals(order='postorder'):
	n.tip_count = sum([c.tip_count for c in n])
	for c in n:
		c.up=n

vic_lineage = T.common_ancestor([leaves['B/Shanghai/1/77'],
								 leaves['B/Victoria/02/1987']])

Phylo.write(vic_lineage, out_vic, 'newick')

n = vic_lineage.up
n_prev = vic_lineage
while n.tip_count == vic_lineage.tip_count:
	n_prev = n
	n=n.up

n.clades = [c for c in n if c!=n_prev]

Phylo.write(T, out_yam, 'newick')
