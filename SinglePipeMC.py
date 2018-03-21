import pylab as pp
import numpy as np
import os

from Transient_Calcs import *
Directory = 'Projects/SinglePipeMC/'

FileName = 'SinglePipeTest.inp'
Net = Import_EPANet_Geom(Directory+FileName)
Iterations = 100
Output = np.zeros((23184,Iterations))
for i in range(Iterations):
	#Net.geom_Plot(plot_Node_Names = True)
	Net.Import_EPANet_Results()
	Net.set_fixed_friction(np.random.normal(0.1,0.001))
	Wavespeed = 1159.2
	Nodes = 20
	Net.Constant_Wavespeed(Wavespeed)
	dt = 50/(Wavespeed*Nodes)
	Net.MOC_Initialisation(dt)


	#Net.Assign_Emmiters('Emitters.csv')

	#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
	Net.Control_Input(Directory+'Demands.csv')
	Net.MOC_Run(50)
	Output[:,i] = Net.nodes[1].TranH
	

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


