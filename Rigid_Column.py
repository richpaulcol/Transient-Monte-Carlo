import numpy as np
import pylab as pp
import wntr as wntr
import networkx as nx
from Transient_Calcs import *

def create_Incidence_Matrix(Net):
	NodesNum = len(Net.nodes)
	LinksNum = len(Net.links)
	Net.Type0 = []
	Net.Type1 = []
	for i in Net.nodes:
		if i.type == 'Reservoir':
			Net.Type0.append(i.number)
		else:
			Net.Type1.append(i.number)

	Net.A = np.zeros((NodesNum,LinksNum))
	Net.A0 = np.zeros((len(Net.Type0),LinksNum))
	Net.A1 = np.zeros((len(Net.Type1),LinksNum))
	for i in Net.links:
		LinkNumber = i.number
		LinkStart = i.node1.number
		LinkEnd = i.node2.number
		Net.A[LinkStart,LinkNumber]=1
		Net.A[LinkEnd,LinkNumber] = -1

	count = 0
	for i in Net.Type0:
		Net.A0[count,:]=Net.A[i,:]
		count+=1

	count = 0
	for i in Net.Type1:
		Net.A1[count,:] = Net.A[i,:]
		count+=1

def create_Solution_Vectors(Net):
	NodesNum = len(Net.nodes)
	LinksNum = len(Net.links)
	Net.Q = np.zeros((LinksNum,1))
	Net.q0 = np.zeros((len(Net.Type0),1))
	Net.q1 = np.zeros((len(Net.Type1),1))
	Net.H0 = np.zeros((len(Net.Type0),1))
	Net.H1 = np.zeros((len(Net.Type1),1))
	
def create_Jacobian_submatrix(Net):
	LinksNum = len(Net.links)
	Net.D = np.zeros((LinksNum,LinksNum))
	if Net.Friction_Units == 'H-W':
		n = 1.852
		for i in range(LinksNum):
			#if Net.links[i].type == 'Pipe'
			K = 10.67*Net.links[i].length / (Net.links[i].diameter**4.87 * Net.links[i].roughness)
			Net.D[i,i] = n * K * abs(Net.Q[i])**(n-1)
	

Directory = 'Projects/Real_Networks/'
FileName = 'ExampleNetwork1.inp'

Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()
create_Incidence_Matrix(Net)
create_Solution_Vectors(Net)
create_Jacobian_submatrix(Net)	
	
