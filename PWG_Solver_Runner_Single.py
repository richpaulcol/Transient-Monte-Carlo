import pylab as pp
import numpy as np
import os

from Transient_Calcs import *

def demand_generator(filename,node,maxTime,dt,originalFlow,maxFlow,startTime,endTime):
	"""
	This function will create a simple transient demand control file for a single node
	It accepts as arguments:
	'filename' this is the required filename of the output file.
	'node' this is the node name that you want to specify a control on
	'maxTime' this is the simulation maximum time in seconds
	'dt' the simulation dt
	'originalFlow' the original flow at that node as specified in the SS simulation
	'maxFlow' the flow you require after the flow change
	
	"""
	times = np.arange(0,maxTime,dt)
	demands = np.ones(times.shape)*originalFlow
	demands[int(startTime/dt):int(endTime/dt)] = maxFlow
	np.savetxt(filename,np.hstack((np.array(node),demands)),delimiter = ',')

Directory = 'Projects/SinglePipe/'

FileName = 'SinglePipe.inp'
Net = Import_EPANet_Geom(Directory+FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(300)
Net.MOC_Initialisation(0.05)

#Net.Assign_Emmiters('Emitters.csv')

demand_generator(Directory+'single_driver.csv',2,150,0.05,1,2,2,150)
#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
Net.Control_Input(Directory+'single_driver.csv')
Net.MOC_Run(150)
Net.transient_Node_Plot()


