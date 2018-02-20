import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = '5Pipes.inp'
Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(300)
Net.MOC_Initialisation(0.05)

#Net.Assign_Emmiters('Emitters.csv')

#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
Net.Control_Input('Kalman_Trial.csv')
Net.MOC_Run(150)
#Net.MOC_Run(86350)
#Net.geom_Plot(plot_Node_Names = True)
#Net.transient_Node_Plot(['6','10','13','16','21','24','31'])

Net.transient_Node_Plot(['1','2','3','4','5','6'])

#for Node in Net.nodes:
#	np.save('Deterministic_Solutions_as_inputs/MeasureData'+str(Node.Name)+'.npy',Node.TranH)
	
#for Node in Net.nodes:
#	pp.scatter([int(Node.Name)], [np.mean(Node.TranH)])
#PE = np.zeros(9999)#999)
#KE = np.zeros(9999)
#for pipe in Net.pipes:
#	PE += np.array(pipe.PE)
#	KE += np.array(pipe.KE)


