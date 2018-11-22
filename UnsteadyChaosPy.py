import numpy as np 
import pylab as pp
import chaospy as cp
Directory = 'Projects/UncertainShear/'
#QData =np.load(Directory+'PCQData.npy')
#HData = np.load(Directory+'PCHData.npy')
#ShearData = np.load(Directory+'PCShearData.npy')
TurbData = np.load(Directory+'PCTurbData.npy')
Samples = np.load(Directory+'PCSamples.npy')

Alpha = cp.Normal(1,0.1)
Beta = cp.Normal(0.1,0.01)
F = cp.Normal(0.02,0.001)
K = cp.Uniform(0.01,0.05)
Distributions = cp.J(Alpha,Beta,F,K)


maxT = 60.0
dT = 0.002
time = np.arange(0,maxT+dT,dT)

Order = 3
NoSamples = 40
Output = TurbData[:NoSamples,:,-1]
polynomial_expansion = cp.orth_ttr(Order, Distributions)
foo_approx = cp.fit_regression(polynomial_expansion, Samples[:,:NoSamples], Output[:,-1])
expected = cp.E(foo_approx, Distributions)
deviation = cp.Std(foo_approx, Distributions)
COV = cp.Cov(foo_approx, Distributions)
#Perc = cp.Perc(foo_approx,[5,95],Distributions)

f,axs = pp.subplots(figsize=(9, 6),nrows = 1,ncols = 1,sharex=True)
axs.plot(time[1:],expected,'k')
axs.fill_between(time[1:],expected+deviation,expected-deviation,color='k',alpha=0.25)
#axs.fill_between(time[1:],Perc[0],Perc[1],colour = 'k',alpha= 0.25)

pp.show()



