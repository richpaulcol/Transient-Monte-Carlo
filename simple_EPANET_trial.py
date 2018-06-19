#import numpy as np
#import pylab as pp
from epanettools import epanet2 as epa
import csv
import math
import numpy as np
#import os

def open_EPANET(filename):
	ret = epa.ENopen(filename,filename[:-3]+'rep',filename[:-3]+'out')

def run_EPANET():
	ret = epa.ENopenH()
	ret = epa.ENinitH(0)
	epa.ENrunH()
	
def argmax(Data):
	Max = max(Data)
	

Directory = 'Projects/SafeStandOff/'		## File location
FileName = 'Net2LPS.inp'			## File name

FileName = str(input("Enter EPANET .inp filename (please include full path i.e. /Directory/Filename.inp):"))

open_EPANET(Directory+FileName)			## Opening the file using the epanettools
run_EPANET()					## initial run of the epanet file

ret,no_nodes=epa.ENgetcount(epa.EN_NODECOUNT)	## Getting the number of nodes in the file
	
Heads = {'Node':'Head'}				##Initialising a dictionary to store the Head data in
Pressures = {'Node':'Pressure'}			##Initialising a dictionary to store the Pressure data in
Demands = {'Node':'Demand'}			##Initialising a dictionary to store the Demand data in

for i in range(1,no_nodes+1):
	ret,index=epa.ENgetnodeid(i)		## Getting the node name
	
	ret,H0 = epa.ENgetnodevalue(i, epa.EN_HEAD )		## getting the nodal head results
	ret,P0 = epa.ENgetnodevalue(i,epa.EN_PRESSURE)		## getting the nodal pressure results
	ret,Demand = epa.ENgetnodevalue(i, epa.EN_DEMAND)	## getting the nodal demand results
	
	Heads[index] = H0					## putting the nodal head results into the heads dictionary
	Pressures[index] = P0					## putting the nodal pressure results into the pressures dictionary
	Demands[index] = Demand					## putting the nodal demand results into the demands dictionary

f = open(FileName[:-3] + 'Initial.csv','wb')
writer = csv.DictWriter(f, Heads.keys())
writer.writeheader()
writer.writerows((Heads,Pressures,Demands))
f.close()


Heads = {'Node':'Head'}
Pressures = {'Node':'Pressure'}
Demands = {'Node':'Demand'}
WSSCd = {'Node':'WSSC'}
NEN3650d = {'Node':'NEN3650_3'}
Orifice = {'Node':'Orifice Size'}

maxD = 100.
Cd = 0.6



for i in range(1,no_nodes+1):
	ret,index=epa.ENgetnodeid(i)
	
	Hs = []
	Ps = []
	Qs = []
	Rbs = []
	dOs = []
	for j in range(100):
		dO = j/maxD
		AO = dO**2*math.pi/4
		ret = epa.ENsetnodevalue(i,epa.EN_EMITTER,1000*Cd*AO*math.sqrt(2*9.81))
		
		
		run_EPANET()
		ret,H = epa.ENgetnodevalue(i, epa.EN_HEAD )		## getting the nodal head results
		ret,P = epa.ENgetnodevalue(i,epa.EN_PRESSURE)		## getting the nodal pressure results
		ret,Demand = epa.ENgetnodevalue(i, epa.EN_DEMAND)	## getting the nodal demand results
		ret,Base = epa.ENgetnodevalue(i,epa.EN_BASEDEMAND)
		
		Q_failure = abs(Demand - Base)
		
		if dO != 0:
			Hs.append(H)
			Qs.append(Q_failure)
			Ps.append(P)
			
			Rbs.append(1.8 *(Q_failure * H/(9.81**0.5 * Cd * dO**(7./2.)))**0.243 * dO)
			dOs.append(dO)
		else:
			WSSCd[index] = (1/1.7)*(Cd**2*H)*(1-math.tan(math.radians(30)))
		#print node.Name,node.P_0,get_epanet_final_demand(node.Name),Rb
		
	
	Heads[index] = Hs[Rbs.index(max(Rbs))]
	Pressures[index] = Ps[Rbs.index(max(Rbs))]
	Demands[index] = Qs[Rbs.index(max(Rbs))]
	NEN3650d[index] = Rbs[Rbs.index(max(Rbs))]
	Orifice[index] = dOs[Rbs.index(max(Rbs))]
	ret = epa.ENsetnodevalue(i,epa.EN_EMITTER,0)
	

f = open(FileName[:-3] + 'Crater_Size.csv','wb')
writer = csv.DictWriter(f, Heads.keys())
writer.writeheader()
writer.writerows((Heads,Pressures,Demands,WSSCd,NEN3650d,Orifice))
f.close()



