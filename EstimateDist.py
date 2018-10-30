import numpy as np
import pylab as pp
from scipy.stats import skewnorm
from scipy.optimize import fmin

Directory = 'Projects/Pressure_Surges_Models/'
MCMean = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Mean.npy',mmap_mode='r')
MCStd = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Std.npy',mmap_mode='r')
MCP5 = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Percent5.npy',mmap_mode='r')
MCP95 = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Percent95.npy',mmap_mode='r')

Time = 1600

M = MCMean[0,Time]  
S = MCStd[0,Time]**2
P5 = MCP5[0,Time]
P95 = MCP95[0,Time]
 


def alphaFind(alpha):
	xi = (M*np.pi*alpha**2-S**(1./2.)*(np.pi*alpha**2-2*alpha**2+np.pi)**(1./2.)*2**(1./2.)*alpha-2*M*alpha**2+M*np.pi)/(np.pi*alpha**2-2*alpha**2+np.pi)
	omega = 1/(np.pi*alpha**2-2*alpha**2+np.pi)**(1./2.)*S**(1./2.)*(alpha**2+1)**(1./2.)*np.pi**(1./2.)
	Dave = skewnorm(alpha, xi,omega)
	
	return (P5-Dave.ppf([0.05]))**2 + (P95-Dave.ppf([0.95]))**2
alpha = fmin(alphaFind,0)


xi = (M*np.pi*alpha**2-S**(1./2.)*(np.pi*alpha**2-2*alpha**2+np.pi)**(1./2.)*2**(1./2.)*alpha-2*M*alpha**2+M*np.pi)/(np.pi*alpha**2-2*alpha**2+np.pi)
omega = 1/(np.pi*alpha**2-2*alpha**2+np.pi)**(1./2.)*S**(1./2.)*(alpha**2+1)**(1./2.)*np.pi**(1./2.)
Dave = skewnorm(alpha, xi,omega)

mean, var, skew, kurt = Dave.stats(moments='mvsk')

print P5/Dave.ppf([0.05]),P95/Dave.ppf([0.95]),M,mean,S,var

#x = np.linspace(Dave.ppf(0.001),Dave.ppf(0.999), 100)
#pp.plot(x, Dave.pdf(x),'r-', lw=5, alpha=0.6, label='skewnorm pdf')
pp.figure()
pp.hist(Dave.rvs(10000),100,density=True)
pp.xlabel('Pressure Head (m)')
pp.ylabel('Probability')
pp.title('Time: '+str(Time*0.01) + ' (s)')

pp.show()
pp.savefig(Directory+'MCDist'+str(Time*0.01)+'.png')

Samples = Dave.rvs(10**6)
####  Simulated Convergence
x = np.arange(2,6.1,0.1)
x2 = 10**x
O = []
for i in x2:
	O.append(np.std(Samples[:int(i)])**2)

pp.semilogx(x2,O)
