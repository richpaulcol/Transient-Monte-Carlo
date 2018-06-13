import pylab as pp
import numpy as np
import pykalman as pk
from Transient_Calcs import pyhyd

#def PODDS_Advect(TMass,IncomingTMass,Q,tau_a,tau_s,alpha,beta):
def PODDS_Advect(State):
	tau_a = State[0]
	tau_s = State[1]
	beta = 2.1

	#tau_a = pyhyd.shear_stress(D, QQ, PipeRoughness, T=10.0, den=1000.0, force_turb_flow=True)
	#if tau_s<tau_a:
	#print tau_a-tau_s
	tau_s = tau_s + dt * beta*(tau_a - tau_s)
	
	return np.hstack((tau_a,tau_s))


def PODDS_Measure(State):
	return State[0]

### Physical Properties
D = 0.1				#Pipe Diameter (m)
PipeArea = np.pi*D**2 / 4	#Pipe CSA (m**2)
Qmax = 0.008			#Max Flow Rate (l/s)
PipeLength = 10		#Pipe Length (m)
PipeRoughness = 0.0001		#Pipe Roughness (m)
maxT = 300			#Simulation Time (s)
#tau_s = 0.00			#Initial Condition shear (N/m^2)
alpha = 5.			#Initial Guess of amount fo material per m ^2
beta = 0.6			#Initial Guess of mobilisation rate (1/s)
Sol1 = alpha/beta
Sol2 = alpha/beta**0.5

### Simulation Properties
n = 2.				#Resolution Tweak (-)
dt = 1.0/n			#Simulation Time Step (s)
Umax = Qmax / PipeArea		#Max Velocity (m/s)
dx = dt * Umax			#Spatial Discretisation stepsize (m)
x = np.arange(0,PipeLength,dx)	#Spatial Discretisation 
TMass = np.ones(x.shape)*1e-4	#Initial Turbidity (kg)
times = np.arange(0,maxT+dt,dt) #Simulation Time Steps

Q = np.ones(times.size)*0.001
#Q[int(50/dt):] = 0.002
#Q[int(100/dt):] = 0.005
#Q[int(200/dt):] = 0.006
#Q[int(200/dt):int(220/dt)] = 0.005 + (np.arange(int(200/dt),int(220/dt)) - (200/dt)) * (0.006-0.005) /(20/dt)
#Q[int(250/dt):] = 0.007

#Q[int(50/dt):int(250/dt)] = 0.001 + (np.arange(int(50/dt),int(250/dt)) - (50/dt)) * (0.006-0.001) /(200/dt)
tau_s = pyhyd.shear_stress(D, Q[0], PipeRoughness, T=10.0, den=1000.0, force_turb_flow=True)



State = np.hstack((Q[0],0))
States = np.zeros((State.size,times.size))
States[:,0] = State
Covariances = np.zeros((State.size,State.size,times.size))

Std = np.ones(States.shape[0])
Std[0] = 1e-1
Std[1] = 1e-1

CC  = np.array([[1e-1,0],[0,1e-1]])

Covariances[:,:,0] =  np.diag(Std)



Obs_Covariance = np.array([1e-4])

Trans_Variance = np.ones(States.shape[0])
Trans_Variance[0] = 1e-5
Trans_Variance[1] = 1e-4			#Transition Uncertainty in Beta
Trans_Covariance = CC#+ np.diag(Trans_Variance) #np.ones((States.shape[0],States.shape[0]))/100000. 




AUSKF = pk.AdditiveUnscentedKalmanFilter(PODDS_Advect,PODDS_Measure,transition_covariance =Trans_Covariance,observation_covariance = Obs_Covariance, initial_state_mean=States[:,0],initial_state_covariance = Covariances[:,:,0])


QMeasured = Q+ np.random.normal(0,1e-4,times.size)


MeasuringTimeInterval =10 #(s)
MeasurementInterval = int(MeasuringTimeInterval/dt)



for t in range(times.size)[1:]:
	
	#if t%MeasurementInterval == 0:
	States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],Q[t])
		#print t*dt
#		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],np.hstack((TMeasured[t],T2Measured[t],QMeasured[t])).T)
	#else:
	#States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1])
	
pp.figure()
pp.subplot(211)
pp.plot(times,Q)
pp.plot(times,States[0,:])
pp.fill_between(times,States[0,:]-Covariances[0,0,:],States[0,:]+Covariances[0,0,:],alpha=0.5)

pp.subplot(212)
pp.plot(times,States[1,:])
pp.fill_between(times,States[1,:]-Covariances[1,1,:],States[1,:]+Covariances[1,1,:],alpha=0.5)

pp.show()

