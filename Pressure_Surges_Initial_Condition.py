import pylab as pp
import numpy as np
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
	pp.figure(figsize=(10,6))
	pp.plot(times,np.mean(Data,axis=1),'k')
	pp.fill_between(times,np.mean(Data,axis=1)+np.std(Data,axis=1),np.mean(Data,axis=1)-np.std(Data,axis=1),alpha = 0.5,color='k',linestyle = '--')
#	pp.fill_between(times,np.percentile(Data,25,axis=1),np.percentile(Data,75,axis=1),alpha = 0.5,color='k',linestyle = '--')
	pp.fill_between(times,np.percentile(Data,5,axis=1),np.percentile(Data,95,axis=1),alpha = 0.25,color='k',linestyle = ':')
	pp.ylabel('Head (m)')
	pp.xlabel('Time (s)')
	pp.xlim(times[0],times[-1])
	pp.tight_layout()
	pp.show()

Directory = 'Projects/Pressure_Surges_Models/'
FileName = '5Pipes.inp'

Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()

Iterations = 1000000
#Inputs
UpstreamHead = np.random.normal(20,0.1,Iterations)
DownstreamFlow = np.random.normal(0.001,0.0001,Iterations)
Demand = np.random.exponential(0.0005,Iterations)

Frictions = np.random.normal(1,0.4,(5,Iterations))

ret,reservoir = epa.ENgetnodeindex(str(Net.nodes[-1].Name))
ret,demand = epa.ENgetnodeindex(str(Net.nodes[1].Name))
ret,downstream = epa.ENgetnodeindex(str(Net.nodes[-2].Name))

Outputs = np.zeros((Iterations,11))

NodeList = [Net.nodes[-1]]+Net.nodes[:-1]
for i in range(Iterations):
	##	Setting the upstream head
	ret = epa.ENsetnodevalue(reservoir,epa.EN_ELEVATION,UpstreamHead[i])
	##	Setting the mid point demand
	ret = epa.ENsetnodevalue(demand,epa.EN_BASEDEMAND,Demand[i]*1000.)
	##	Setting the downstream initial flow
	ret = epa.ENsetnodevalue(downstream,epa.EN_BASEDEMAND,DownstreamFlow[i]*1000.)
	#print 'Dave',ret
	
	for j in range(len(Net.pipes)):
		ret,index = epa.ENgetlinkindex(str(Net.pipes[j].Name))
		ret,epa.ENsetlinkvalue(index,epa.EN_ROUGHNESS,Frictions[j,i])
	
	Net.run_epanet_file()
	
	OutputState = []
	
	
	for pipe in Net.pipes:
		ret,index = epa.ENgetlinkindex(str(pipe.Name))
		ret,Flow =  epa.ENgetlinkvalue(index,epa.EN_FLOW)
#		print pipe.Name,Flow
		OutputState.append(Flow)
		
	for node in NodeList:
		ret,index = epa.ENgetnodeindex(str(node.Name))
		ret,Head =  epa.ENgetnodevalue(index,epa.EN_HEAD)
		#print node.Name,Head
		OutputState.append(Head)
	
	Outputs[i,:] = OutputState
	
np.save(Directory+'InitialConditions.npy',Outputs)
