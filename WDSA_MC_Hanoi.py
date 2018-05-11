import pylab as pp
import numpy as np
import chaospy as cp
import os
from Transient_Calcs import *

def normal_dist(mean,variance):
	stddev = np.sqrt(variance)
	x = np.linspace(mean-3*stddev,mean+3*stddev,100)
	y = 1./(np.sqrt(2 * np.pi * variance)) * np.exp(- (x-mean)**2 / (2*variance))
	return x,y

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

def output_plotter(times,Data):
	pp.figure()
	pp.plot(times,np.mean(Data,axis=1),'b')
	pp.fill_between(times,np.mean(Data,axis=1)+2*np.std(Data,axis=1),np.mean(Data,axis=1)-2*np.std(Data,axis=1),alpha = 0.5,color='b')
	pp.show()


Directory = 'Projects/WDSA_Models/'
FileName = '5Pipes.inp'
FileName = 'hanoi3.inp'
#FileName = 'Net3LPSnoControl.inp'
Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()


distribution = cp.J(cp.Normal(1.5e-5,1e-6))

Wavespeed = 1000.
dt =0.01	

Iterations = 1000
samples = distribution.sample(Iterations)
#Output = np.zeros((4637,Iterations))
maxTime = 50

demand_generator(Directory+'Hanoi_driving_transient.csv',13,maxTime,dt,1,2,0.5,2)
Transient_Times = np.arange(0,maxTime,dt)

Output = np.zeros((Iterations,len(Net.nodes),Transient_Times.size))

for i in range(Iterations):
	print i
#	#Net.geom_Plot(plot_Node_Names = True)
	Net.alter_epanet_friction(samples[i])		## This function assigns a new friction value to the pipes on each iteration
	Net.run_epanet_file()				## This function runs the SS EPAnet model	
	Net.read_results_from_epanet()			## This function reads the SS EPAnet data into the transinet model
	Net.Assign_Emmiters_All()			## This function assigns all the demands as emiters in the transient model (i.e. pressure based demands) with a CdA value chosen to match the EPAnet demand.
	Net.Constant_Wavespeed(Wavespeed)		## This function sets a constant wavespeed for all pipes in the model

	Net.MOC_Initialisation(dt)			## This function initialises the transient model and discretises the pipe lengths with a constant dt

	Net.Control_Input(Directory+'Hanoi_driving_transient.csv')		## This function applies the given demand based control to the model using the specified file.
	Net.MOC_Run(maxTime)				## This is the function that actually runs the MOC simulation
	count = 0					## The next few lines of code reads the output of the nodal pressures and feeds them into the Output array specified earlier.
	for node in Net.nodes:
		Output[i,count,:] = node.TranH
		count+=1
#	Output[:,i] = Net.nodes[1].TranH
##	Output[0,i] = Net.nodes[1].H_0
#	
Net.close_epanet_file()					## Ensuring that any EPAnet file is closed at the end of the run (basic housekeeping)

np.save(Directory + 'Hanoi_varying_friction.npy',Output)			## Saving the output data



