import numpy as np
import pylab as pp
from scipy import signal
maxT = 60
dT = 0.002
dX = 1
D = 0.05
A = np.pi*D**2/4
time = np.arange(0,maxT+dT,dT)
NoNodes = 101

Directory = 'Projects/UncertainShear/'
LinearHeadFast = np.load(Directory+'LinearHeadFast.npy')
LinearFlowFast = np.load(Directory+'LinearFlowFast.npy')
LinearShearFast = np.load(Directory+'LinearShearFast.npy')

LinearHeadSlow = np.load(Directory+'LinearHeadSlow.npy')
LinearFlowSlow = np.load(Directory+'LinearFlowSlow.npy')
LinearShearSlow =np.load(Directory+'LinearShearSlow.npy')

USHeadFast = np.load(Directory+'USHeadFast.npy')
USFlowFast = np.load(Directory+'USFlowFast.npy')
USShearFast = np.load(Directory+'USShearFast.npy')

USHeadSlow = np.load(Directory+'USHeadSlow.npy')
USFlowSlow = np.load(Directory+'USFlowSlow.npy')
USShearSlow = np.load(Directory+'USShearSlow.npy')

#f,axs = pp.subplots(3,1,sharex=True)
#axs[0].plot(time,LinearHeadFast[:,NoNodes/2],label = 'Linear Fast')
#axs[0].plot(time,LinearHeadSlow[:,NoNodes/2],label = 'Linear Slow')
#axs[0].plot(time,USHeadFast[:,NoNodes/2],label = 'US Fast')
#axs[0].plot(time,USHeadSlow[:,NoNodes/2],label = 'US Slow')

#axs[1].plot(time,LinearFlowFast[:,NoNodes/2],label = 'Linear Fast')
#axs[1].plot(time,LinearFlowSlow[:,NoNodes/2],label = 'Linear Slow')
#axs[1].plot(time,USFlowFast[:,NoNodes/2],label = 'US Fast')
#axs[1].plot(time,USFlowSlow[:,NoNodes/2],label = 'US Slow')

#axs[2].plot(time,LinearShearFast[:,NoNodes/2],label = 'Linear Fast')
#axs[2].plot(time,LinearShearSlow[:,NoNodes/2],label = 'Linear Slow')
#axs[2].plot(time,USShearFast[:,NoNodes/2],label = 'US Fast')
#axs[2].plot(time,USShearSlow[:,NoNodes/2],label = 'US Slow')


#axs[0].set_ylabel('Head (m)')
#axs[1].set_ylabel('Flow (m3/s)')
#axs[2].set_ylabel('Wall Shear (Pa)')

#axs[0].legend()
#pp.show()

### Advect
alpha = 5.			#Initial Guess of amount fo material per m ^2
beta = .1
Shear = LinearShearFast
Flow  = LinearFlowFast

dX2 = 0.01
TMass = np.zeros(int(100/dX2))
tau_s = Shear[0]
tau_s_list = []
TMass_list = []
PODDSMass_list = []
dN_list = []

for t in range(1,len(time)):
	tau_s = tau_s + dT * beta*(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	dN = alpha * beta *(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	dN2 = signal.resample(dN, int(100/dX2))
	
	PODDSMass = np.pi*D*dX2*dN2
	lam = ( (signal.resample(Flow[t], int(100/dX2))/A)*dT / dX2)
	TMass[1:] = lam[:-1]*TMass[:-1] + (1-lam[1:])*TMass[1:] + PODDSMass[1:]
	TMass[0] = 0
	
	PODDSMass_list.append(PODDSMass)
	dN_list.append(dN)
	tau_s_list.append(tau_s)
	TMass_list.append(TMass)
tau_s_A = np.array(tau_s_list)
TMassA = np.array(TMass_list)
	
