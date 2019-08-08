'''
Script comparing different approximations of the susceptibility to the full effect
of the product of all previous infections.
'''

import numpy as np
import matplotlib.pyplot as plt

xmax=3
xoa = np.linspace(0,xmax,101)
alpha = 1.0
R0 = 2.5
eps = 1.0 #alpha/np.log(R0)*0.7

def K(x, alpha):
	'''
	cross-immunity induced by an infection at distance x where
		alpha is the protection against homotypic infection
	'''
	return 1-alpha*np.exp(-x)


def full_product(x, mmax=10, alpha=1.0, eps=1.0):
	'''
	product of susceptibilities of periodic infections.
		mmax is the maximum number of terms
		alpha is the degree of protection against homotypic infection
		eps is the reinfection interval in units of the cross-immunity scale
	'''
	S = np.ones_like(x)
	for m in range(mmax):
		S*=K(x+m*eps, alpha)
	return S


def R_R(x, alpha):
	'''
	Rouzine and Rozhnova's approximation of dropping all but the last term
	'''
	return K(x,alpha)

def linearized(x, alpha=1.0, eps=1.0):
	'''
	full linearization of the log of the product
	'''
	return np.exp(-alpha*np.exp(-x)*np.exp(eps)/(np.exp(eps)-1))

def first_p_linearized(x, alpha=1.0, eps=1.0):
	'''
	keep the first term, linear and resum all other terms
	'''
	return K(x,alpha)*np.exp(-alpha/(np.exp(eps)-1)*np.exp(-x))


## generate panel A: compare different approxiations as a function of antigenic distance
fs=16
fig, axs = plt.subplots(1,2,figsize=(12,6))
axs[0].plot(xoa, -np.log(full_product(xoa, alpha=alpha, eps=eps)), label=r'$S(x)$',lw=3)
axs[0].plot(xoa, -np.log(linearized(xoa, alpha=alpha, eps=eps)), label='weak inhibition',lw=2)
#axs[0].plot(xoa, -np.log(first_p_linearized(xoa, alpha=alpha, eps=eps)), label='linear2')
axs[0].plot(xoa, -np.log(R_R(xoa,alpha=alpha)), label='most recent',lw=2)
# axs[0].plot([0,xmax], [np.log(R0), np.log(R0)])
# axs[0].plot([eps,eps], [0,3])
axs[0].legend(fontsize=0.8*fs)
axs[0].set_yscale('log')
axs[0].set_xlabel(r'antigenic distance of last infection $x/d$', fontsize=fs)
axs[0].set_ylabel(r'neg. log susceptibility $-\log S(x)$', fontsize=fs)
axs[0].tick_params(labelsize=0.8*fs)
axs[0].text(-0.15, 0.9, '(A)', transform=axs[0].transAxes, fontsize=1.5*fs)

## integrate periodic infection model over the population distribution:
# when parameterized by the last infection, the population is uniform with density 1/d eps up to d eps
# the overall recovered density is flat 1/d eps out to infinity

def R_R_total(alpha, eps):
	'''
	here, we have a constant density in the interval [0,eps] split into
	100 intervals for numeric integration
	'''
	return np.mean(R_R(np.linspace(0,eps,101), alpha))


def full_model(alpha, eps):
	'''
	here, we have a constant density in the interval [0,eps] split into
	100 intervals for numeric integration
	'''
	return np.mean(full_product(np.linspace(0,eps,101), alpha=alpha, eps=eps))


def MF(alpha, eps):
	'''
	In the mean field model, we simply integrate over all the density of recovered,
	which is uniform and extends out to -infinity
	'''
	return np.exp(-alpha/eps)

## generate panel B of the figure
eps_range = np.linspace(0.3, 1, 20)
for ls,alpha in zip([ '--', '-'], [0.75, 1.0]):
	axs[1].plot(eps_range, [full_model(alpha, y) for y in eps_range], label=r'full model, $\alpha='+str(alpha)+'$', ls=ls, c='C0', lw=3)
	axs[1].plot(eps_range, MF(alpha, eps_range), label=r'mean-field, $\alpha='+str(alpha)+'$', ls=ls, c='C1', lw=2)
	axs[1].plot(eps_range, [R_R_total(alpha, y) for y in eps_range], label=r'most recent,  $\alpha='+str(alpha)+'$', ls=ls, c='C2', lw=2)

plt.legend(fontsize=0.8*fs, loc=2)
axs[1].set_ylabel('population average susceptibility', fontsize=fs)
axs[1].set_xlabel(r'antigenic distance between reinfections $\epsilon$', fontsize=fs)
axs[1].tick_params(labelsize=0.8*fs)
axs[1].set_ylim(0,0.8)
axs[1].text(-0.15, 0.9, '(B)', transform=axs[1].transAxes, fontsize=1.5*fs)
plt.tight_layout()
plt.savefig('susceptibility.png')

# generate additional figures for different alpha
for alpha in [0.5, 0.8, 1.0]:
	plt.figure()
	plt.title(f"alpha={alpha}")
	plt.plot(eps_range, [full_model(alpha, y) for y in eps_range], label='full_model', lw=3)
	plt.plot(eps_range, [R_R_total(alpha, y) for y in eps_range], label='R&R')
	plt.plot(eps_range, MF(alpha, eps_range), label='MF')

	plt.legend()
	plt.savefig(f"model_comparison_alpha_{alpha}.png")
