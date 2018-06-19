import pylab as pp
import numpy as np
from Transient_Calcs import *
#import operator


def normal_dist(mean,variance):
	stddev = np.sqrt(variance)
	x = np.linspace(mean-3*stddev,mean+3*stddev,100)
	y = 1./(np.sqrt(2 * np.pi * variance)) * np.exp(- (x-mean)**2 / (2*variance))
	return x,y


def set_epanet_emitter(node,value):
	ret,index = epa.ENgetnodeindex(str(node))
	ret = epa.ENsetnodevalue(index,epa.EN_EMITTER,value)

def get_epanet_final_demand(node):
	ret,index = epa.ENgetnodeindex(str(node))
	ret,value = epa.ENgetnodevalue(index,epa.EN_DEMAND)
	return value
def steady_pressure_bar_chart():
	Pressures = {}
	for node in Net.nodes:
		Pressures[str(node.Name)] = node.H_0
	pp.figure()	
	pp.bar(range(len(Pressures)), list(Pressures.values()), align='center')
	pp.xticks(range(len(Pressures)), list(Pressures.keys()))
	pp.show()

def find_local_high_points():
	HighPoints=[]
	for node in Net.nodes:
		ThisNodeElevation = node.zPos
		OtherNodeElevations = []
		for pipe in node.pipesIn:
			OtherNodeElevations.append(pipe.node1.zPos)
		for pipe in node.pipesOut:
			OtherNodeElevations.append(pipe.node2.zPos)
		if (ThisNodeElevation>=np.array(OtherNodeElevations)).all():
			HighPoints.append(node.Name)

	return HighPoints
			
		

Cd = 1.

Directory = 'Projects/SafeStandOff/'
FileName = 'Net2LPS.inp'
#FileName = 'hanoi3.inp'
#FileName = 'Net3LPSnoControl.inp'
Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()

Heads = {}
Flows = {}
PitSize = {}
Orifice = {}
PreFlows ={}
Net.run_epanet_file()
Net.read_results_from_epanet()
for pipe in Net.pipes:
	PreFlows[pipe.Name] = abs(pipe.Q_0 * 0.06309 * 1000.)

WSSCd = {}
for node in Net.nodes:
	Hs = []
	Qs = []
	Rbs = []
	dOs = []
	for dO in np.linspace(0,0.5,100):
		AO = dO**2*np.pi/4
		set_epanet_emitter(node.Name,1000*Cd*AO*np.sqrt(2*9.81))
		Net.run_epanet_file()
		Net.read_results_from_epanet()
		H = node.P_0 
		Q = get_epanet_final_demand(node.Name)
		Hs.append(H)
		Qs.append(Q)
		Rbs.append(1.8 *(Q * H/(9.81**0.5 * Cd * dO**(7./2.)))**0.243 * dO)
		dOs.append(dO)
		if dO == 0:
			WSSCd[node.Name] = (1/1.7)*(Cd**2*H)*(1-np.tan(np.radians(30)))
		#print node.Name,node.P_0,get_epanet_final_demand(node.Name),Rb
		
	
	Heads[node.Name] = Hs[np.argmax(Rbs[1:])+1]
	Flows[node.Name] = Qs[np.argmax(Rbs[1:])+1]
	PitSize[node.Name] = Rbs[np.argmax(Rbs[1:])+1]
	Orifice[node.Name] = dOs[np.argmax(Rbs[1:])+1]
	set_epanet_emitter(node.Name,0)	 
	

#Net.close_epanet_file()

#diameters = []
#for pipe in Net.pipes:
#	diameters.append(pipe.diameter*1000)

pp.figure()
	
pp.bar(range(len(Heads)), list(Heads.values()), align='center',alpha = 0.5,color='r')
pp.xticks(range(len(Heads)), list(Heads.keys()))
pp.ylabel('Head (m)')
pp.ylim(0,80)
#pp.xlabel('Node')
#pp.show()

#pp.figure()	
pp.twinx()
pp.bar(range(len(Flows)), list(Flows.values()), align='center',alpha = 0.5)
pp.xticks(range(len(Flows)), list(Flows.keys()))
pp.ylabel('Leak Flow (l/s)')
pp.ylim(400,0)
pp.xlabel('Node')
pp.show()

pp.figure()
pp.bar(range(len(PitSize)), list(PitSize.values()), align='center',alpha = 0.5)
pp.xticks(range(len(PitSize)), list(PitSize.keys()))
pp.bar(range(len(WSSCd)), list(WSSCd.values()), align='center',alpha = 0.5)
pp.xticks(range(len(WSSCd)), list(WSSCd.keys()))
#pp.ylim(40,0)
pp.ylabel('Pit Size (m)')
pp.xlabel('Pipe')
pp.show()
#Highpoints = find_local_high_points()

