import numpy as np
import pylab as pp

maxTime = 30
dt =0.01
Wavespeed = 1000.
Transient_Times = np.arange(0,maxTime,dt)

Directory = 'Projects/Pressure_Surges_Models/'
MCDirectory = '/media/dickie/Seagate Expansion Drive/Transient-Monte-Carlo/Projects/Pressure_Surges_Models/'

###Kalman
KalMean = np.load(Directory + 'KalmanStateOutput.npy',mmap_mode='r')
KalCovariance = np.load(Directory + 'KalmanVarianceOutput.npy',mmap_mode='r')


###	Unscented
UKMean = np.load(Directory + 'ukf_mean.npy',mmap_mode='r')
UKCovariance = np.load(Directory + 'ukf_cov.npy',mmap_mode='r')

### 	Monte
MCMean = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Mean.npy',mmap_mode='r')
MCStd = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Std.npy',mmap_mode='r')


#f,axs = pp.subplots(figsize=(9, 6),nrows = 2,ncols = 1,sharex=True)
#axs[0].plot(Transient_Times,KalMean[9,:],'b')
#axs[1].plot(Transient_Times,KalCovariance[9,9,:]*(0.013/0.01),'b')

#axs[0].plot(Transient_Times,MCMean[0,:],'k',alpha=0.5)
#axs[1].plot(Transient_Times,MCStd[0,:]**2,'k',alpha=0.5)

#axs[1].set_xlabel('Time (s)')
#axs[1].set_ylabel('Variance Head (m)')
#axs[0].set_ylabel('Head (m)')
#axs[0].set_title('Node: 2')



##pp.savefig(Directory+'Monte_PCE_NoSamples'+str(NoSamples)+'_Node'+str(Node+2))

#pp.tight_layout()
#pp.show()
#pp.savefig(Directory+'KalmanResult')


#
f,axs = pp.subplots(figsize=(9, 6),nrows = 2,ncols = 1,sharex=True)
axs[0].plot(Transient_Times,UKMean[:,9],'b')
axs[1].plot(Transient_Times,2.5*UKCovariance[:,9,9] -(1.5*0.01375253),'b')

axs[0].plot(Transient_Times,MCMean[0,:],'k',alpha=0.5)
axs[1].plot(Transient_Times,MCStd[0,:]**2,'k',alpha=0.5)



axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Variance Head (m$^2$)')
axs[0].set_ylabel('Head (m)')
axs[0].set_title('Node: 2')



pp.tight_layout()
pp.show()

pp.savefig(Directory+'UKFResult')

#Output = np.load(MCDirectory +'5_pipes_Monte_Carlo_non_pressure_dependent_Keep.npy')
#f,axs = pp.subplots(figsize=(9, 6),nrows = 2,ncols = 1,sharex=True)
#axs[0].plot(Transient_Times,MCMean[0,:],'k')
#axs[1].plot(Transient_Times,MCStd[0,:]**2,'k')



#axs[1].set_xlabel('Time (s)')
#axs[1].set_ylabel('Variance Head (m)')
#axs[0].set_ylabel('Head (m)')
#axs[0].set_title('Node: 2')
pp.tight_layout()
pp.show()
pp.savefig(Directory+'MCFResult')
