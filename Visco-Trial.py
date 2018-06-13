import numpy as np
import pylab as pp
D = 0.1
alpha = 0.1
gamma = 0.1
h = 0.01

tMax = 1000
N = 10
dt = 0.1


time = np.arange(0,tMax)
H = np.zeros(time.shape)
H[1:] = 1
H[100:] = 10.

tau = np.ones(N)*np.exp(np.arange(N))
J = np.ones(N)*30.

y = np.zeros((N,tMax))
strain = np.zeros((tMax))
for t in time[:-1]:
	for i in range(N):
		y[i,t+1] = np.exp(-dt / tau[i]) * y[i,t] + J[i]/tau[i] * np.exp(-dt/tau[i]) *(H[t+1] - H[t])
	strain[t+1] = alpha * D / (2 * h) * gamma *sum(y[:,t+1])

pp.figure()
pp.plot(time,H)
pp.plot(time,strain)

pp.show()
