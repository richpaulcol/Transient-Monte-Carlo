import pylab as pp
import numpy as np
import os

from Transient_Calcs import *

Directory = 'Projects/State_Parameter_Estimation_Real/'

FileName = 'Real_Branched.inp'
Net = Import_EPANet_Geom(Directory+FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(1000)
Net.MOC_Initialisation(0.00001)
Net.Assign_Emmiters_All()
Net.MOC_Run(1)
#Net.Assign_Emmiters('Emitters.csv')

#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
#Net.Control_Input(Directory+'Kalman_Trial.csv')
#Net.MOC_Run(150)
#Net.MOC_Run(86350)
#Net.geom_Plot(plot_Node_Names = True)
#Net.transient_Node_Plot(['6','10','13','16','21','24','31'])

#Net.transient_Node_Plot(['1','2','3','4','5','6'])

#for Node in Net.nodes:
#	np.save(Directory +'MeasureData'+str(Node.Name)+'.npy',Node.TranH)
	
#for Node in Net.nodes:
#	pp.scatter([int(Node.Name)], [np.mean(Node.TranH)])
#PE = np.zeros(9999)#999)
#KE = np.zeros(9999)
#for pipe in Net.pipes:
#	PE += np.array(pipe.PE)
#	KE += np.array(pipe.KE)


