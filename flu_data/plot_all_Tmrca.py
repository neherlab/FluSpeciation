import matplotlib.pyplot as plt
import numpy as np

lineages = [('H3N2', 'A/H3N2'), ('H1N1pdm', 'A/H1N1pdm(2009)'), ('H1N1', 'A/H1N1(1918)'),
			('Bvic','B/Vic'),
			('Byam', 'B/Yam')]
fs=20
cutoff = 1970
ymax=15
dTmrca = {}
mTmrca = {}
Tmrca = {}
for l,ln in lineages:
	try:
		tmrca_traj = np.loadtxt('%s_tmrca_trajectory.dat'%l)
		tmrca_traj = tmrca_traj[tmrca_traj[:,0]>cutoff]
		if tmrca_traj.shape[0]%2==1:
			tmrca_traj = tmrca_traj[1:]
		dTmrca[ln] = np.maximum(0,np.diff(tmrca_traj[:,1])[::2])
		mTmrca[ln] = np.sum(0.5*(np.maximum(0,tmrca_traj[1::2,1]) + np.maximum(0,tmrca_traj[0::2,1]))*(tmrca_traj[1::2,0] - tmrca_traj[0::2,0]))
		mTmrca[ln]/= tmrca_traj[-1,0] - tmrca_traj[0,0]
		Tmrca[ln] = tmrca_traj
	except Exception as e:
		print("can't load",ln, e)

fig, axs = plt.subplots(5,1, figsize=(12,16), sharex=True)
for ax,(l, ln) in zip(axs, lineages):
	if ln in Tmrca:
		tmrca_traj = Tmrca[ln]
		ax.plot(tmrca_traj[:,0], tmrca_traj[:,1], lw=3)
		ax.text(cutoff+1.5, ymax*0.5,
				ln + "\n" + r"$\langle T_{MRCA}\rangle=" + "%1.1f$"%mTmrca[ln]
				   + "\n" + r"$\langle \delta T_{MRCA}\rangle="+"%1.1f$"%dTmrca[ln].mean(), fontsize=fs)
	ax.set_ylabel(r'$T_{MRCA}$', fontsize=fs)
	ax.set_ylim(0,15)
	# ax.legend(loc=2, fontsize=0.8*fs)
	ax.tick_params(labelsize=0.8*fs)

plt.xlabel('year', fontsize=fs)
plt.savefig('all_Tmrca.pdf')
