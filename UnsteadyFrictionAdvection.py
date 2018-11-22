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
LinearShearFast = abs(np.load(Directory+'LinearShearFast.npy'))

LinearHeadSlow = np.load(Directory+'LinearHeadSlow.npy')
LinearFlowSlow = np.load(Directory+'LinearFlowSlow.npy')
LinearShearSlow = abs(np.load(Directory+'LinearShearSlow.npy'))

USHeadFast = np.load(Directory+'USHeadFast.npy')
USFlowFast = np.load(Directory+'USFlowFast.npy')
USShearFast = abs(np.load(Directory+'USShearFast.npy'))

USHeadSlow = np.load(Directory+'USHeadSlow.npy')
USFlowSlow = np.load(Directory+'USFlowSlow.npy')
USShearSlow = abs(np.load(Directory+'USShearSlow.npy'))

#f,axs = pp.subplots(3,1,sharex=True,figsize = [9,9])
#axs[0].plot(time,LinearHeadFast[:,NoNodes/2],label = 'Linear Fast')
#axs[0].plot(time,LinearHeadSlow[:,NoNodes/2],label = 'Linear Slow')
##axs[0].plot(time,USHeadFast[:,NoNodes/2],label = 'US Fast')
##axs[0].plot(time,USHeadSlow[:,NoNodes/2],label = 'US Slow')

#axs[1].plot(time,LinearFlowFast[:,NoNodes/2],label = 'Linear Fast')
#axs[1].plot(time,LinearFlowSlow[:,NoNodes/2],label = 'Linear Slow')
##axs[1].plot(time,USFlowFast[:,NoNodes/2],label = 'US Fast')
##axs[1].plot(time,USFlowSlow[:,NoNodes/2],label = 'US Slow')

#axs[2].plot(time,LinearShearFast[:,NoNodes/2],label = 'Linear Fast')
#axs[2].plot(time,LinearShearSlow[:,NoNodes/2],label = 'Linear Slow')
##axs[2].plot(time,USShearFast[:,NoNodes/2],label = 'US Fast')
##axs[2].plot(time,USShearSlow[:,NoNodes/2],label = 'US Slow')


#axs[0].set_ylabel('Head (m)')
#axs[1].set_ylabel('Flow (m3/s)')
#axs[2].set_ylabel('Wall Shear (Pa)')
#axs[2].set_xlabel('Time (s)')
#axs[0].legend()
#pp.tight_layout()
#pp.show()
##pp.xlim(0,10)
#pp.savefig(Directory+'LinearComp.png')



#sfsfsfsadfasdf

### Advect
alpha = 1.			#Initial Guess of amount fo material per m ^2
beta = .1
Shear = USShearFast
Flow  = USFlowFast

dX2 = 0.01
TMass = np.zeros(NoNodes)
tau_s = Shear[0]
tau_s_list = []
TMass_list = []
PODDSMass_list = []
dN_list = []

for t in range(1,len(time)):
	tau_s = tau_s + dT * beta*(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	dN = alpha * beta *(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	#dN2 = signal.resample(dN, int(100/dX2))
	
	PODDSMass = np.pi*D*dX*dN
	
	TMassNew = np.zeros(TMass.shape)
#	if t%1 ==0:
	lam = ( Flow[t]/A)*dT / dX
	TMassNew[1:] = lam[:-1]*TMass[:-1] + (1-lam[1:])*TMass[1:] + PODDSMass[1:]
#	else: 
#		TMassNew = TMass+PODDSMass
	
	
	PODDSMass_list.append(PODDSMass)
	dN_list.append(dN)
	tau_s_list.append(tau_s)
	TMass_list.append(TMassNew)
	TMass = TMassNew
	
	
SSTime = np.arange(0,10*maxT+dT,dT)
for tt in range(SSTime.size):
	tau_s = tau_s + dT * beta*(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	dN = alpha * beta *(Shear[t] - tau_s)*((Shear[t] - tau_s)>0)
	PODDSMass = np.pi*D*dX*dN
	
	TMassNew = np.zeros(TMass.shape)
#	if t%1 ==0:
	lam = ( Flow[t]/A)*dT / dX
	TMassNew[1:] = lam[:-1]*TMass[:-1] + (1-lam[1:])*TMass[1:] + PODDSMass[1:]
#	else: 
#		TMassNew = TMass+PODDSMass
	
	
	PODDSMass_list.append(PODDSMass)
	dN_list.append(dN)
	tau_s_list.append(tau_s)
	TMass_list.append(TMassNew)
	TMass = TMassNew

tau_s_A = np.array(tau_s_list)
TMassA = np.array(TMass_list)

TotalTime = np.arange(0,TMassA.shape[0]*dT,dT)
pp.plot(TotalTime,TMassA[:,-1],label = 'USFast')
