import numpy as np
import pylab as pp

dt = 0.01
Steps = 1000

Customers = 100

k = 1		### The base poisson coefficient
sigma = 1	### The standard deviation of the delta demand
Demands = np.zeros((Customers,Steps))
for t in range(1,Steps):
	for i in range(Customers):
		deltaDemand = np.random.poisson(k*dt)*np.random.normal(-1.75,sigma)
#		print deltaDemand
		if Demands[i,t-1] + deltaDemand > 0:
			Demands[i,t] = Demands[i,t-1] + deltaDemand
		else:
			Demands[i,t] = 0
			
pp.plot(Demands.sum(axis=0))
pp.show()
