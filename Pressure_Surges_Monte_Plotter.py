import numpy as np
import pylab as pp
Directory = 'Projects/Pressure_Surges_Models/'

Mean = np.load(Directory + '5_pipes_Monte_Carlo_Mean.npy')
Std = np.load(Directory + '5_pipes_Monte_Carlo_Std.npy')
Percent5 = np.load(Directory + '5_pipes_Monte_Carlo_Percent5.npy')
Percent95 = np.load(Directory + '5_pipes_Monte_Carlo_Percent95.npy')


maxTime = 30
dt =0.01
Wavespeed = 1000.
Transient_Times = np.arange(0,maxTime,dt)

Node = 2


pp.figure()
pp.plot(Transient_Times,Mean[Node,:],'b')
#pp.fill_between(Transient_Times, Mean[Node,:]+Std[Node,:],Mean[Node,:]-Std[Node,:],color='b',alpha = 0.5)
pp.fill_between(Transient_Times, Percent5[Node,:],Percent95[Node,:],color='b',alpha = 0.25)



pp.show()
