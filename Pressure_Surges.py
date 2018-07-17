import pylab as pp
import numpy as np
import chaospy as cp
import os
from Transient_Calcs import *
import pykalman as pk
#import seaborn as sns

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

Directory = 'Projects/Pressure_Surges_Models/'
FileName = '5Pipes.inp'
maxTime = 300
dt =0.01
Wavespeed = 1000.
Transient_Times = np.arange(0,maxTime,dt)

#####
# 	The input distributions
#	Upstream Head =  N(20,0.1)  [m]
#	Downstream Flow = N(0.001,0.0001)   [m3/s]
#	Demand at node3 = Exp(0.0005)
#	Roughness = N(1,0.4)

Upstream = cp.Normal(20,0.1)
Downstream = cp.Normal(0.001,0.0001)
Demand = cp.Exponential(0.0005)
Roughness = cp.J(cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4))

Distributions = cp.J(Upstream,Downstream,Demand,Roughness)
#samples = Distributions.sample(100000,'S')
#np.save(Directory+'Transient_MC_samples.npy',samples)


#
#####		Monte Carlo Method
#
#

#samples = np.load(Directory+'Transient_MC_samples.npy')

#Net = Import_EPANet_Geom(Directory+FileName)
#Net.open_epanet_file()


#Iterations = 1
#Output = np.zeros((Iterations,len(Net.nodes),Transient_Times.size))

#ret,reservoir = epa.ENgetnodeindex(str(Net.nodes[-1].Name))
#ret,demand = epa.ENgetnodeindex(str(Net.nodes[1].Name))
#ret,downstream = epa.ENgetnodeindex(str(Net.nodes[-2].Name))

#for i in range(Iterations):
#	print i
#	ret = epa.ENsetnodevalue(reservoir,epa.EN_ELEVATION,samples[0,i])
#	ret = epa.ENsetnodevalue(downstream,epa.EN_BASEDEMAND,samples[1,i]*1000.)
#	ret = epa.ENsetnodevalue(demand,epa.EN_BASEDEMAND,samples[2,i]*1000.)
#	for j in range(len(Net.pipes)):
#		ret,index = epa.ENgetlinkindex(str(Net.pipes[j].Name))
#		ret,epa.ENsetlinkvalue(index,epa.EN_ROUGHNESS,samples[3+j,i])
#		
#	Net.run_epanet_file()				## This function runs the SS EPAnet model	
#	Net.read_results_from_epanet()			## This function reads the SS EPAnet data into the transinet model
#	demand_generator(Directory+'5_pipes_driving_transient.csv',6,maxTime,dt,samples[1,i]*1000.,samples[1,i]*1000.+0.5,2,30)
#	#Net.Assign_Emmiters_All()
#	Net.Constant_Wavespeed(Wavespeed)
#	Net.MOC_Initialisation(dt)
#	Net.Control_Input(Directory+'5_pipes_driving_transient.csv')
#	Net.MOC_Run(maxTime)
#	
#	count = 0					## The next few lines of code reads the output of the nodal pressures and feeds them into the Output array specified earlier.
#	for node in Net.nodes:
#		Output[i,count,:] = node.TranH
#		count+=1
		
#np.save(Directory + '5_pipes_Monte_Carlo_non_pressure_dependent.npy',Output)

######		Linear Transient
##
##	Calculating the mean inputs

Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()

####	Setting up the initial condition
ret,reservoir = epa.ENgetnodeindex(str(Net.nodes[-1].Name))
ret,demand = epa.ENgetnodeindex(str(Net.nodes[1].Name))
Net.nodes[1].demand = 0.0005
ret,downstream = epa.ENgetnodeindex(str(Net.nodes[-2].Name))
Net.nodes[-2].demand = 0.001
ret = epa.ENsetnodevalue(reservoir,epa.EN_ELEVATION,20)
ret = epa.ENsetnodevalue(downstream,epa.EN_BASEDEMAND,0.001*1000.)
ret = epa.ENsetnodevalue(demand,epa.EN_BASEDEMAND,0.0005*1000.)
for j in range(len(Net.pipes)):
	ret,index = epa.ENgetlinkindex(str(Net.pipes[j].Name))
	ret,epa.ENsetlinkvalue(index,epa.EN_ROUGHNESS,1)

#####	Running the initial steady state
Net.run_epanet_file()				## This function runs the SS EPAnet model	
Net.read_results_from_epanet()
Net.Constant_Wavespeed(Wavespeed)
Net.Initialise_Linear_Kalman(dt)
Net.Assign_Emmiters_All()
Net.dx = Wavespeed*dt
####	Initialising the Covariance Matrices
COV = np.load(Directory+'InitialCOV.npy')
Var = np.diag(COV)

Net.P_Matrix = Net.P_Matrix.todense()
Net.Q_Matrix = Net.Q_Matrix.todense()
#for i in range(Net.pipes_State_Index.size):
#	for j in range(Net.pipes_State_Index.size):
#		k = Net.pipes_State_Index.astype('int')[i]
#		l = Net.pipes_State_Index.astype('int')[j]
#		Net.P_Matrix[Net.CPs+i,Net.CPs+j] = COV[k,l]
#		Net.P_Matrix[i,j] = COV[k+len(Net.pipes),l+len(Net.pipes)]
		
		
#Net.P_Matrix[:Net.CPs,:Net.CPs] = (np.ones((50,50))*0.15)
#Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs] = (np.ones((50,50))*0.01)
#Net.P_Matrix[Net.CPs:Net.CPs+20,Net.CPs:Net.CPs+20] = (np.ones((20,20))*0.26)
Net.P_Matrix[Net.CPs-1,Net.CPs-1] = 0.1  		#Variance in the upstream Head BC#
Net.P_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.0001**2	#Variance in the downstream flow#
Net.P_Matrix[2*Net.CPs+1,2*Net.CPs+1] = 0.0005**2	#Variance in the node 3 demand

Net.Q_Matrix[Net.CPs-1,Net.CPs-1] = 0.1  		#Variance in the upstream Head BC#
Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.0001**2	#Variance in the downstream flow#
Net.Q_Matrix[2*Net.CPs+1,2*Net.CPs+1] = 0.0005**2	#Variance in the node 3 demand
		
#		Net.Q_Matrix[Net.CPs+i,Net.CPs+i] = 8 * Net.dx* Net.pipes[k].Q_0**2 / (9.81*Net.pipes[k].diameter * np.pi**2)*1.7337e-06
#Net.initial_BC_Uncertainty_Head(0.1)
kf = pk.KalmanFilter()

kf.transition_matrices = Net.A_Matrix.todense()

kf.transition_covariance = Net.Q_Matrix

kf.initial_state_mean = Net.X_Vector
kf.initial_state_covariance = Net.P_Matrix
Simon = kf.transition_matrices
Dave = Net.regenerateA()
StateOutput = np.zeros((Net.X_Vector.size,Transient_Times.size))
VarianceOutput = np.zeros((Net.X_Vector.size,Net.X_Vector.size,Transient_Times.size))
for i in range(0,Transient_Times.size):
	if i == int(2/dt):
		Net.X_Vector[2*Net.CPs+Net.nodes[-2].number]+=0.0005
		
	Net.X_Vector,Net.P_Matrix = kf.filter_update(Net.X_Vector,Net.P_Matrix)
	StateOutput[:,i] = Net.X_Vector
	VarianceOutput[:,:,i] = Net.P_Matrix 
#	Dave = Net.regenerateA()
#	kf.transition_matrices = Net.regenerateA()
#	

np.save('StateOutput.npy',StateOutput)
np.save('VarianceOutput.npy',VarianceOuput)

import Pressure_Surges_Monte_Plotter
f,axs = pp.subplots(nrows = 2,ncols = 1)
axs[0].plot(Transient_Times,StateOutput[Net.nodal_CPs['3'],:])
axs[1].plot(Transient_Times,VarianceOutput[Net.nodal_CPs['3'],Net.nodal_CPs['3'],:])

#f,axs = pp.subplots(nrows = 4,ncols = 2)
#axs[0,0].imshow(Simon[:50,:50])
#axs[0,1].imshow(kf.transition_matrices[:50,:50])

#axs[1,0].imshow(Simon[50:,50:])
#axs[1,1].imshow(kf.transition_matrices[50:,50:])

#axs[2,0].imshow(Simon[:50,50:])
#axs[2,1].imshow(kf.transition_matrices[:50,50:])

#axs[3,0].imshow(Simon[50:,:50])
#axs[3,1].imshow(kf.transition_matrices[50:,:50])




pp.show()
