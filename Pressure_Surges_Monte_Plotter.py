import numpy as np
import pylab as pp
Directory = 'Projects/Pressure_Surges_Models/'

Mean = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Mean.npy')
Std = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Std.npy')
Percent5 = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Percent5.npy')
Percent95 = np.load(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent_Percent95.npy')


maxTime = 30
dt =0.01
Wavespeed = 1000.
Transient_Times = np.arange(0,maxTime,dt)

Node = 1

f,axs = pp.subplots(figsize=(9, 6),nrows = 2,ncols = 1,sharex=True)
axs[0].plot(Transient_Times,Mean[Node,:],'b')
axs[1].plot(Transient_Times,Std[Node,:]**2,'b')

axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Variance Head (m)')
axs[0].set_ylabel('Head (m)')
axs[0].set_title('Node 2')
pp.tight_layout()
#pp.figure()
#pp.plot(Transient_Times,Mean[Node,:],'b')
#pp.fill_between(Transient_Times, Mean[Node,:]+Std[Node,:],Mean[Node,:]-Std[Node,:],color='b',alpha = 0.5)
#pp.fill_between(Transient_Times, Percent5[Node,:],Percent95[Node,:],color='b',alpha = 0.25)



pp.show()
