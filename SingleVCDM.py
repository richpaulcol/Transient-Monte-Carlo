import pylab as pp
import numpy as np
import pykalman as pk
from Transient_Calcs import pyhyd

#def PODDS_Advect(TMass,IncomingTMass,Q,tau_a,tau_s,alpha,beta):
def PODDS_Advect(State):
	TMass = State[:-5]
	IncomingTMass = State[-5]
	tau_a = State[-4]
	tau_s = State[-3],
	alpha = State[-2] 
	beta = State[-1]

#	if beta <0:
#		beta = 0
#	if alpha <0:
#		alpha = 0
#	if Q < 0:
#		Q = 0
#	if IncomingTMass <0:
#		IncomingTMass = 0

	#tau_a = pyhyd.shear_stress(D, Q, PipeRoughness, T=10.0, den=1000.0, force_turb_flow=False)
	#if tau_s<tau_a:
	tau_s = tau_s + dt * beta*(tau_a - tau_s)
	dN = alpha * beta *(tau_a - tau_s)	
	#else:
	#tau_s = tau_s + dt * beta*(tau_a - tau_s) / 1000.
	#dN = alpha * beta *(tau_a - tau_s)	


	PODDSMass = np.pi*D*dx*dN
	lam = ( (Q[t]/PipeArea)*dt / dx)
	#if Q>0:
	TMass[1:] = lam*TMass[:-1] + (1-lam)*TMass[1:] + PODDSMass
	TMass[0] = IncomingTMass + PODDSMass
	#else:
		
	#print dN,tau_a-tau_s
	return np.hstack((TMass,IncomingTMass,tau_a,tau_s,alpha,beta))

#def PODDS_Measure(TMass,IncomingTMass,Q,tau_a,tau_s,alpha,beta):
def PODDS_Measure(State):
	TMass = State[:-5]
	#IncomingTMass = State[-6]
	tau_a = State[-4]
#	tau_a = State[-4]
#	tau_s = State[-3],
#	alpha = State[-2] 
#	beta = State[-1]
#	return np.hstack((TMass[-1],TMass[x.size/2],tau_a))
	return np.hstack((TMass[x.size/2],tau_a))

### Physical Properties
D = 0.1				#Pipe Diameter (m)
PipeArea = np.pi*D**2 / 4	#Pipe CSA (m**2)
Qmax = 0.01			#Max Flow Rate (l/s)
PipeLength = 10		#Pipe Length (m)
PipeRoughness = 0.0001		#Pipe Roughness (m)
maxT = 500			#Simulation Time (s)
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

#Driving Functions (
IncomingTMass = np.ones(times.size)*0.05#np.random.uniform(0,0.1,times.size)
IncomingTMass[int(150/dt):] = 0.25
IncomingTMass[int(350/dt):] = 0.25 + (np.arange(301)) * (0.0-0.25) /(301)
IncomingTMass += np.random.normal(0,IncomingTMass/10.,times.size)

Q = np.ones(times.size)*0.001
Q[int(50/dt):] = 0.002
Q[int(100/dt):] = 0.005
Q[int(200/dt):] = 0.006
Q[int(200/dt):int(220/dt)] = 0.005 + (np.arange(int(200/dt),int(220/dt)) - (200/dt)) * (0.006-0.005) /(20/dt)
Q[int(250/dt):] = 0.007
Q[int(400/dt):] = 0.008

#Q[int(50/dt):int(250/dt)] = 0.001 + (np.arange(int(50/dt),int(250/dt)) - (50/dt)) * (0.006-0.001) /(200/dt)

tau_a = pyhyd.shear_stress(D, Q, PipeRoughness, T=10.0, den=1000.0, force_turb_flow=False)
tau_s = pyhyd.shear_stress(D, Q[0], PipeRoughness, T=10.0, den=1000.0, force_turb_flow=False)

TMass0 = TMass
tau_s0 = tau_s
T = [0]
T2 = [0]
tau_s0s = [tau_s]
for t in range(times.size)[1:]:	
	State = np.hstack((TMass0,IncomingTMass[t],tau_a[t],tau_s0,alpha,beta))
	State = PODDS_Advect(State)
	tau_s0 = State[-3]
	TMass0 = State[:-5]
	T.append(State[-6])
	T2.append(State[x.size/2])
	tau_s0s.append(State[-3])

T = np.array(T)
T2 = np.array(T2)


#Output Functions

alpha = 5.0			#Initial Guess of amount fo material per m ^2
beta = 0.8
State = np.hstack((TMass,IncomingTMass[0],tau_a[0],tau_s,alpha,beta))
States = np.zeros((State.size,times.size))
States[:,0] = State
Covariances = np.zeros((State.size,State.size,times.size))

Std = np.ones(States.shape[0])
Std[:-5] = 1e-4	#Initial Uncertainty in TMass along the pipe length
Std[-5] = 1e-4			#Initial Uncertainty in incomingTMass
Std[-4] = 1e-1 		#Initial Uncertainty in tau_a
Std[-3] = 1e-3			#Initial Uncertainty in tau_s
Std[-2] = 1e-5			#Initial Uncertainty in Alpha
Std[-1] = 1e-1			#Initial Uncertainty in Beta

#Std = np.random.uniform(0,0.1,States.shape[0])

Covariances[:,:,0] =  np.diag(Std)#+np.diag(Std[:-1]/10.,1)+np.diag(Std[:-1]/10.,-1)np.ones((States.shape[0],States.shape[0]))/100000. 
#Obs_Covariance = np.array([[1e-8]])
Obs_Covariance = np.array([[1e-8,0],[0,1e-6]])
#Obs_Covariance = np.array([[1e-2,0,0],[0,1e-2,0],[0,0,1e-8]])

Trans_Variance = np.ones(States.shape[0])
Trans_Variance[:-5] = (dt**2)*1e-2	
Trans_Variance[-5] = (dt**2)*1e-2		#Transition Uncertainty in incomingTMass
Trans_Variance[-4] = (dt**2)*1e-1			#Transition Uncertainty in tau_a
Trans_Variance[-3] = (dt**2)*1e-1			#Transition Uncertainty in tau_s
Trans_Variance[-2] = (dt**2)*1e-10			#Transition Uncertainty in Alpha
Trans_Variance[-1] = (dt**2)*1e-3			#Transition Uncertainty in Beta
Trans_Covariance = + np.diag(Trans_Variance) #np.ones((States.shape[0],States.shape[0]))/100000. 

AUSKF = pk.AdditiveUnscentedKalmanFilter(PODDS_Advect,PODDS_Measure,transition_covariance =Trans_Covariance,observation_covariance = Obs_Covariance, initial_state_mean=States[:,0],initial_state_covariance = Covariances[:,:,0])

#State,Covariance = AUSKF.filter_update(State,Covariance,np.array([0.1,0.1]))
#T[::2] = np.ma.masked
#T[::3] = np.ma.masked
Mask = np.ones(times.size)

#TMeasured = np.ma.masked_array(T+ np.random.normal(0,1e-2,times.size),TMask)
#QMeasured = np.ma.masked_array(Q+ np.random.normal(0,5e-4,times.size),QMask)


TMeasured = T+ np.random.normal(0,1e-2,times.size)
T2Measured = T2 + np.random.normal(0,1e-2,times.size)
tau_Measured = tau_a+ np.random.normal(0,1e-2,times.size)
MeasuringTimeInterval =3 #(s)
MeasurementInterval = int(MeasuringTimeInterval/dt)



for t in range(times.size)[1:]:
	
	if t%MeasurementInterval == 0:
		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],np.hstack((T2Measured[t],tau_Measured[t])))
		#print t*dt
		#States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],np.hstack((TMeasured[t],T2Measured[t],tau_Measured[t])).T)
	else:
		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1])
	
#	States[States[:,t-1] < 0,t-1] = 0




pp.figure()
pp.subplot(411)
pp.plot(times,T,'k')
pp.plot(times,T2,'g')
pp.scatter(times[::MeasurementInterval],TMeasured[::MeasurementInterval],c = 'k')
pp.scatter(times[::MeasurementInterval],T2Measured[::MeasurementInterval],c = 'g')
pp.plot(times,States[-6,:],'r')
pp.fill_between(times,States[-6,:]+Covariances[-6,-6,:],States[-6,:]-Covariances[-6,-6,:],color = 'r',alpha = 0.5)
pp.plot(times,States[x.size/2,:],'y')
pp.fill_between(times,States[x.size/2,:]+Covariances[x.size/2,x.size/2,:],States[x.size/2,:]-Covariances[x.size/2,x.size/2,:],color = 'y',alpha = 0.5)
#pp.plot(times,States[0,:])
#pp.plot(times,States[20,:])
pp.subplot(412)
pp.plot(times,tau_a,'k')
pp.scatter(times[::MeasurementInterval],tau_Measured[::MeasurementInterval],c='k')
pp.plot(times,States[-4,:],'r')
pp.fill_between(times,States[-4,:]+Covariances[-4,-4,:],States[-4,:]-Covariances[-4,-4,:],color = 'r',alpha = 0.5)
pp.subplot(413)
pp.plot(times,tau_s0s,'k')
pp.plot(times,States[-3,:],'r')
pp.fill_between(times,States[-3,:]+Covariances[-3,-3,:],States[-3,:]-Covariances[-3,-3,:],color = 'r',alpha = 0.5)
pp.subplot(414)
pp.plot(times,IncomingTMass,'k')
pp.plot(times,States[-5,:],'r')
pp.fill_between(times,States[-5,:]+Covariances[-5,-5,:],States[-5,:]-Covariances[-5,-5,:],color = 'r',alpha = 0.5)
#pp.subplot(515)
#pp.plot(times,States[0,:])
#pp.plot(times,States[4,:])
#pp.plot(times,States[8,:])
#pp.plot(times,States[12,:])
pp.show()

np.save('States.npy',States)
np.save('Covariances.npy',Covariances)
np.save('Times.npy',times)
np.save('Incoming.npy',IncomingTMass)
np.save('tau_s.npy',tau_s0s)
np.save('Q.npy',tau_a)
np.save('Turbidity.npy',T)

print Sol1,Sol2,2*Sol1,2*Sol2
print States[-1,-1],States[-2,-1],States[-2,-1]/States[-1,-1],States[-2,-1]*States[-1,-1]**2
	
