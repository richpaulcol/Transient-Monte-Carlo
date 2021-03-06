import numpy as np
from filterpy.kalman import unscented_transform, MerweScaledSigmaPoints, JulierSigmaPoints, JulierSigmaPoints
import scipy.stats as stats
import pylab as pp

def function(State):
	return np.sin(State)#+np.random.normal(0,1e-4)
#	return State**2

State = np.array([0.1,0.5,0.75,1.0])


ukf_mean = State
ukf_cov = np.identity(State.size)*0.1

#sigmas_f = np.empty((Iterations,sigmas.shape[0],sigmas.shape[1]))
points = MerweScaledSigmaPoints(n=State.size, alpha=0.3, beta=2., kappa=0)
#points = JulierSigmaPoints(n=State.size,kappa=3-State.size)
#points =  JulierSigmaPoints(n=State.size)
means = []
covs = []
states = []
for i in range(100):
	sigmas = points.sigma_points(ukf_mean,ukf_cov)
	sigmas_f = np.empty(sigmas.shape)
	for i in range(sigmas.shape[0]):
		sigmas_f[i] = function(sigmas[i])
	ukf_mean, ukf_cov = unscented_transform(sigmas_f, points.Wm,points.Wc)	
	#ukf_cov+=1e-9
	State = function(State)
	means.append(ukf_mean)
	covs.append(ukf_cov)
	states.append(State)
#	ukf_cov = nearestPD(ukf_cov)

means = np.array(means)
covs = np.array(covs)[:,0,0]
states = np.array(states)


