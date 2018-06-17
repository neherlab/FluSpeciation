import numpy as np
import matplotlib.pyplot as plt
from augur.utils import read_node_data
import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description = "Analyze TMRCA.")
parser.add_argument("--tree", help="tree file")
parser.add_argument("--node_data", help="node_data file")
parser.add_argument("--titers", help="titer_model file")
parser.add_argument("--output", help="output prefix")
args = parser.parse_args()

T = Phylo.read(args.tree, 'newick')
of = [args.node_data, args.titers] if args.titers else [args.node_data]
node_data = read_node_data(of)

T.root.up = None
for n in T.find_clades(order='postorder'):
	n.numdate = node_data["nodes"][n.name]["numdate"]
	if args.titers:
		n.cTiter = node_data["nodes"][n.name]["cTiter"]
		n.dTiter = node_data["nodes"][n.name]["dTiter"]
	if n.is_terminal():
		n.ntips = 1
		n.tree_length = n.branch_length
		if args.titers:
			n.antigenic_length = n.dTiter
	else:
		n.ntips = np.sum([c.ntips for c in n])
		n.tree_length = n.branch_length + np.sum([c.tree_length for c in n])
		if args.titers:
			n.antigenic_length = n.dTiter + np.sum([c.antigenic_length for c in n])
		for c in n:
			c.up=n

to_prune = []
cutoff = 0.03
for n in T.find_clades(order='postorder'):
	if n.tree_length>cutoff and n.ntips<5 and n.up!=T.root:
		print(n.tree_length, n.ntips)
		to_prune.append(n)

for n in to_prune:
	if n.up:
		n.up.clades = [c for c in n.up if c!=n]
		for t in n.get_terminals():
			print("pruning leaf:",t.name)

for n in T.find_clades(order='postorder'):
	if n.is_terminal():
		n.youngest = n.numdate
		n.ntips = 1
		n.tree_length = n.branch_length
		if args.titers:
			n.max_cTiter = n.cTiter
			n.antigenic_length = n.dTiter
	else:
		n.youngest = np.max([c.youngest for c in n])
		n.ntips = np.sum([c.ntips for c in n])
		n.tree_length = n.branch_length + np.sum([c.tree_length for c in n])
		if args.titers:
			n.max_cTiter = np.max([c.max_cTiter for c in n])
			n.antigenic_length = n.dTiter + np.sum([c.antigenic_length for c in n])

nodes_to_bridge = [n for n in T.find_clades() if len(n.clades)==1]
for n in nodes_to_bridge:
	c = n.clades[0]
	p = n.up
	c.up = p
	p.clades = [c] + [k for k in p if n!=k]
	p.clades.sort(key = lambda x:x.ntips)
	c.branch_length += n.branch_length
	if args.titers:
		c.dTiter += n.dTiter


tmrca_traj = []
antigenic_traj = []
antigenic_div_traj = []
n=T.root
n_root = n
loss_date = sorted([c.youngest for c in n])[-2]
#loss_cTiter = sorted([c for c in n], key=lambda x:x.youngest)[-2].max_cTiter
if args.titers:
	loss_cTiter = sorted([c.max_cTiter for c in n])[-2]
while not n.is_terminal():
	if n.up is None:
		tmrca_traj.append((n.numdate, 0))
		antigenic_traj.append((n.numdate, 0))
		antigenic_div_traj.append((n.numdate, 0))
	else:
		new_loss_date = sorted([c.youngest for c in n])[-2]
		if args.titers:
			#new_loss_cTiter = sorted([c for c in n], key=lambda x:x.youngest)[-2].max_cTiter
			new_loss_cTiter = sorted([c.max_cTiter for c in n])[-2]
			antigenic_traj.append((n.numdate, n.cTiter - n_root.cTiter))
			if new_loss_date>loss_date:
				antigenic_div_traj.append((loss_date, sorted([c.max_cTiter for c in n])[-2] - n_root.cTiter))
				antigenic_div_traj.append((loss_date, sorted([c.max_cTiter for c in n])[-2] - n.cTiter))
				# if loss_cTiter<n.cTiter:
				# 	import pdb; pdb.set_trace()
				# antigenic_traj.append((loss_date, loss_cTiter - n_root.cTiter))
				# antigenic_traj.append((loss_date, loss_cTiter - n.cTiter))
				# loss_cTiter = new_loss_cTiter

		if new_loss_date>loss_date:
			tmrca_traj.append((loss_date, loss_date - n_root.numdate))
			tmrca_traj.append((loss_date, loss_date - n.numdate))
			print(loss_date, n.numdate, n.youngest, n.numdate, n_root.numdate)
			print(n.numdate - n_root.numdate)

			loss_date = new_loss_date
			# if (n.numdate<n_root.numdate):
			# 	import ipdb; ipdb.set_trace()
			n_root = n
	next_trunk_i = np.argmax([c.youngest for c in n])
	n = n.clades[next_trunk_i]
tmrca_traj.append((loss_date, loss_date - n_root.numdate))
tmrca_traj = np.array(tmrca_traj)

if args.titers:
	#antigenic_traj.append((n.numdate, sorted([c.max_cTiter for c in n])[-2] - n.cTiter))
	#antigenic_traj.append((loss_date, loss_cTiter - n_root.cTiter))
	antigenic_traj = np.array(antigenic_traj)
	antigenic_div_traj = np.array(antigenic_div_traj)

np.savetxt(args.output + '_tmrca_trajectory.dat', tmrca_traj)
plt.figure(figsize=(12,6))
plt.plot(tmrca_traj[:,0], tmrca_traj[:,1])
plt.ylabel('tree depth')
plt.xlabel('year')
plt.savefig(args.output + '_depth.pdf')

if args.titers:
	# plt.figure()
	np.savetxt(args.output + '_cTiter_trajectory.dat', antigenic_traj)
	plt.plot(antigenic_traj[:,0], antigenic_traj[:,1])
	plt.plot(antigenic_div_traj[:,0], antigenic_div_traj[:,1])
	plt.ylabel('tree depth')
	plt.xlabel('year')
	plt.savefig(args.output + '_titer_depth.pdf')

