import pylab as pp
import numpy as np
import pykalman as pk
from Transient_Calcs import pyhyd
from numpy import linalg as la

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3
    
def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

#def PODDS_Advect(TMass,IncomingTMass,Q,tau_a,tau_s,alpha,beta):
def PODDS_Advect(State):
	TMass = State[:-5]
	IncomingTMass = State[-5]
	tau_a = State[-4]
	tau_s = State[-3]
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
	if tau_s<tau_a:
		dN = dt*alpha * beta *(tau_a - tau_s)	
		#tau_s = tau_s + dt * beta*(tau_a - tau_s)
		#print 'Dave'
		
	else:
		#tau_s = tau_s #+ dt * beta*(tau_a - tau_s) / 1000.
		dN = 0#alpha * beta *(tau_a - tau_s)	
	#print tau_s
	#tau_s = tau_s + dt * beta*(tau_a - tau_s)
	#tau_s -= 0.0001

	PODDSMass = np.pi*D*dx*dN
	lam = ( (Q[t]/PipeArea)*dt / dx)
	#if Q>0:
	TMass[1:] = lam*TMass[:-1] + (1-lam)*TMass[1:] + PODDSMass
	TMass[0] =  PODDSMass #+ IncomingTMass 
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
	
	
Directory = 'Projects/PODDS_Real_Data/'

Data = np.loadtxt(Directory+'SS_Data_500mm_Sept2017.csv',skiprows=1,delimiter=',',dtype='str')
Time = Data[:,0].astype('int')
Flow = Data[:,1].astype('float')
Mask = (Data[:,2]=='')
Data[Mask,2] = 0
Turb = Data[:,2].astype('float')
#Turb = Data[:,2].astype('float')




#Vint = np.cumsum(Velo*dt)

#Steps = []
#for i in range(Vint.size):
#	Steps.append(sum(Vint[i:] - Vint[i]<Length))




### Physical Properties
D = 0.5				#Pipe Diameter (m)
PipeArea = np.pi*D**2 / 4	#Pipe CSA (m**2)
Qmax = np.max(Flow)			#Max Flow Rate (l/s)
PipeLength = 1000		#Pipe Length (m)
PipeRoughness = 0.0001		#Pipe Roughness (m)
StartPoint = 36200
EndPoint = StartPoint + 2000
maxT = 500			#Simulation Time (s)
#tau_s = 0.00			#Initial Condition shear (N/m^2)
alpha = 2.1			#Initial Guess of amount fo material per m ^2
beta = 0.002			#Initial Guess of mobilisation rate (1/s)
Sol1 = alpha/beta
Sol2 = alpha/beta**0.5



Q = Flow[StartPoint:EndPoint+1]
### Simulation Properties
n = 1.				#Resolution Tweak (-)
dt = 300./n			#Simulation Time Step (s)
Umax = np.max(Q) / PipeArea		#Max Velocity (m/s)
dx = dt * Umax			#Spatial Discretisation stepsize (m)
x = np.arange(0,PipeLength,dx)	#Spatial Discretisation 
TMass = np.ones(x.shape)*1e-4	#Initial Turbidity (kg)

times = np.arange(StartPoint*dt,EndPoint*dt+dt,dt) #Simulation Time Steps

#Driving Functions (
IncomingTMass = np.ones(times.size)*0.05#np.random.uniform(0,0.1,times.size)
#IncomingTMass[int(150/dt):] = 1.25
#IncomingTMass[int(350/dt):] = 1.25 + (np.arange(301)) * (0.0-1.25) /(301)
#IncomingTMass += np.random.normal(0,IncomingTMass/10.,times.size)

#Q = np.ones(times.size)*0.001
#Q[int(50/dt):] = 0.002
#Q[int(100/dt):] = 0.005
#Q[int(200/dt):] = 0.006
#Q[int(200/dt):int(220/dt)] = 0.005 + (np.arange(int(200/dt),int(220/dt)) - (200/dt)) * (0.006-0.005) /(20/dt)
#Q[int(250/dt):] = 0.007
#Q[int(400/dt):] = 0.008



#Q[int(50/dt):int(250/dt)] = 0.001 + (np.arange(int(50/dt),int(250/dt)) - (50/dt)) * (0.006-0.001) /(200/dt)

tau_a = pyhyd.shear_stress(D, Q, PipeRoughness, T=10.0, den=1000.0, force_turb_flow=False)
tau_s = tau_a[0]#pyhyd.shear_stress(D, Q[0], PipeRoughness, T=10.0, den=1000.0, force_turb_flow=False)



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
#	T2.append(State[x.size/2])
	tau_s0s.append(State[-3])

T = np.array(T)
#T2 = np.array(T2)


#Output Functions

alpha = 0.1			#Initial Guess of amount fo material per m ^2
beta = 0.002
State = np.hstack((TMass,IncomingTMass[0],tau_a[0],tau_s,alpha,beta))
States = np.zeros((State.size,times.size))
States[:,0] = State
Covariances = np.zeros((State.size,State.size,times.size))

Std = np.ones(States.shape[0])
Std[:-5] = 1e-1			#Initial Uncertainty in TMass along the pipe length
Std[-5] = 1e-2			#Initial Uncertainty in incomingTMass
Std[-4] = 1e-2 			#Initial Uncertainty in tau_a
Std[-3] = 1e-6			#Initial Uncertainty in tau_s
Std[-2] = 1e1			#Initial Uncertainty in Alpha
Std[-1] = 1e1			#Initial Uncertainty in Beta

#Std = np.random.uniform(0,0.1,States.shape[0])

Covariances[:,:,0] =  np.diag(Std)#+np.diag(Std[:-1]/10.,1)+np.diag(Std[:-1]/10.,-1)np.ones((States.shape[0],States.shape[0]))/100000. 
#Obs_Covariance = np.array([[1e-8]])
Obs_Covariance = np.array([[1e-2,0],[0,2e-2]])
#Obs_Covariance = np.array([[1e-2,0,0],[0,1e-2,0],[0,0,1e-3]])

Trans_Variance = np.ones(States.shape[0])
Trans_Variance[:-5] = 1e-2	
Trans_Variance[-5] = 1e-20		#Transition Uncertainty in incomingTMass
Trans_Variance[-4] = 1e1			#Transition Uncertainty in tau_a
Trans_Variance[-3] = 1e-5			#Transition Uncertainty in tau_s
Trans_Variance[-2] = 1e-2			#Transition Uncertainty in Alpha
Trans_Variance[-1] = 1e-2			#Transition Uncertainty in Beta
Trans_Covariance = + np.diag(Trans_Variance) #np.ones((States.shape[0],States.shape[0]))/100000. 

AUSKF = pk.AdditiveUnscentedKalmanFilter(PODDS_Advect,PODDS_Measure,transition_covariance =Trans_Covariance,observation_covariance = Obs_Covariance, initial_state_mean=States[:,0],initial_state_covariance = Covariances[:,:,0])

#State,Covariance = AUSKF.filter_update(State,Covariance,np.array([0.1,0.1]))
#T[::2] = np.ma.masked
#T[::3] = np.ma.masked
Mask = np.ones(times.size)

#TMeasured = np.ma.masked_array(T+ np.random.normal(0,1e-2,times.size),TMask)
#QMeasured = np.ma.masked_array(Q+ np.random.normal(0,5e-4,times.size),QMask)


#TMeasured = T+ np.random.normal(0,1e-2,times.size)
T2Measured = Turb[StartPoint:EndPoint+1]
tau_Measured = tau_a
MeasuringTimeInterval =300 #(s)
MeasurementInterval = int(MeasuringTimeInterval/dt)

Covariances[:,:,0] = nearestPD(Covariances[:,:,0])

for t in range(times.size)[1:]:
	
	if t%MeasurementInterval == 0:
		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],np.hstack((T2Measured[t],tau_Measured[t])))
#		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1],np.hstack((TMeasured[t],T2Measured[t],tau_Measured[t])).T)
	else:
		States[:,t],Covariances[:,:,t] = AUSKF.filter_update(States[:,t-1],Covariances[:,:,t-1])
	
	#rCovariances[:,:,t] = nearestPD(Covariances[:,:,t])
#	States[States[:,t-1] < 0,t-1] = 0




#pp.figure()
#pp.subplot(511)
fig,axs = pp.subplots(3,1,sharex=True)
#pp.plot(times,T,'k')
#pp.plot(times,T2,'g')
#pp.scatter(times[::MeasurementInterval],TMeasured[::MeasurementInterval],c = 'k')
axs[0].scatter(times[::MeasurementInterval],T2Measured[::MeasurementInterval],c = 'g')
axs[0].plot(times,States[-6,:],'r')
axs[0].fill_between(times,States[-6,:]+Covariances[-6,-6,:],States[-6,:]-Covariances[-6,-6,:],color = 'r',alpha = 0.5)
#pp.plot(times,States[x.size/2,:],'y')
#pp.fill_between(times,States[x.size/2,:]+Covariances[x.size/2,x.size/2,:],States[x.size/2,:]-Covariances[x.size/2,x.size/2,:],color = 'y',alpha = 0.5)
#pp.plot(times,States[0,:])
#pp.plot(times,States[20,:])

axs[0].set_ylabel('Turbidity')

#pp.subplot(512)
axs[1].plot(times,tau_a,'k')
axs[1].scatter(times[::MeasurementInterval],tau_Measured[::MeasurementInterval],c='k')
axs[1].plot(times,States[-4,:],'r')
axs[1].fill_between(times,States[-4,:]+Covariances[-4,-4,:],States[-4,:]-Covariances[-4,-4,:],color = 'r',alpha = 0.5)
axs[1].set_ylabel('tau$_a$')
#pp.subplot(513)
#pp.plot(times,tau_s0s,'k')
axs[1].plot(times,States[-3,:],'b')
axs[1].fill_between(times,States[-3,:]+Covariances[-3,-3,:],States[-3,:]-Covariances[-3,-3,:],color = 'r',alpha = 0.5)
axs[1].set_ylabel('tau$_s$')
#pp.subplot(514)
#axs[3].plot(times,IncomingTMass,'k')
#axs[3].plot(times,States[-5,:],'r')
#axs[3].fill_between(times,States[-5,:]+Covariances[-5,-5,:],States[-5,:]-Covariances[-5,-5,:],color = 'r',alpha = 0.5)
#axs[3].set_ylabel('incoming')
#pp.subplot(515)
axs[2].plot(times,States[-1,:])
axs[2].plot(times,States[-2,:])
#pp.plot(times,States[8,:])
#pp.plot(times,States[12,:])
axs[2].set_ylabel(r'$\alpha and \beta$')
pp.show()

#np.save('States.npy',States)
#np.save('Covariances.npy',Covariances)
#np.save('Times.npy',times)
#np.save('Incoming.npy',IncomingTMass)
#np.save('tau_s.npy',tau_s0s)
#np.save('Q.npy',tau_a)
#np.save('Turbidity.npy',T)

#print Sol1,Sol2,2*Sol1,2*Sol2
#print States[-1,-1],States[-2,-1],States[-2,-1]/States[-1,-1],States[-2,-1]*States[-1,-1]**2
	
#fig,axs = pp.subplots(2,1,sharex=True)
#axs[0].plot(times,tau_a,label='tau_a')
##axs[0].plot(times,States[-4,:],label ='tau_a_state')
#axs[0].plot(times,States[-3,:],label ='tau_s_state')
#axs[0].plot(times,tau_s0s,label = 'tau_s')
#axs[0].set_ylim(0,1)

#axs[1].plot(times,T2Measured)
#axs[1].plot(times,States[-6,:])
#axs[1].plot(times,T)


pp.show()
