import pylab as pp
import numpy as np
from PWG_EPA_tools import *
import networkx as nx


FileName = 'Trials/EPANet_Transient_Case_2.inp'
FileName = 'Trials/hanoi3.inp'
#FileName = 'Trials/Dave.inp'
FileName = 'Trials/HampsonModel.inp'
#FileName = 'Trials/Net3.inp'
FileName = 'Trials/YW_B507.inp'
#FileName = 'Trials/Richmond_standard.inp'

Net = Import_EPANet_Geom(FileName)
Net.Import_EPANet_Results()



Net.Constant_Wavespeed(1000)
#Net.MOC_Initialisation(0.000001)
ConnectivityArray = np.zeros((len(Net.nodes),len(Net.nodes)))
TransmissionCoeffArray = np.zeros((len(Net.nodes),len(Net.nodes)))
TimeOfTravelArray = np.zeros((len(Net.nodes),len(Net.nodes)))

for node in Net.nodes:					### Initialising the node information	
	node.TransCoeffs = {}				###	Node Transmission Coefficients Dictionary
	node.pipes = node.pipesIn+node.pipesOut		###	List of all pipes at node
	node.denom = 0					
	for i in range(len(node.pipes)):
		node.denom += node.pipes[i].area / node.pipes[i].c
		
	for i in range(len(node.pipes)):
		node.TransCoeffs[node.pipes[i]] = 2*(node.pipes[i].area / node.pipes[i].c) / node.denom
	#print node.TransCoeffs

G = nx.Graph()
for node in Net.nodes:
	G.add_node(node)
	
for pipe in Net.pipes+Net.valves+Net.pumps:
	pipe.Re = pipe.V_0 * 998 * pipe.diameter / 1e-3
	if pipe.Friction_Units == 'D-W':
		pipe.friction = 0.25 / (np.log10((pipe.roughness/1000.) / (3.7*pipe.diameter) + 5.74/(pipe.Re**0.9))**2) 
		
	elif pipe.Friction_Units == 'H-W':
		pipe.friction = 133.89 / (abs(pipe.V_0)**(4./27.) * pipe.roughness**(50./27.) * pipe.diameter**(1./16.)+10**-10)

	if pipe.Re<2000:
		pipe.friction = 64./(pipe.Re+1)
	
	attenuation = np.exp(-10**2*pipe.friction*pipe.V_0*pipe.length / (2*pipe.c*pipe.diameter))
	#attenuation = pipe.length*np.mean(pipe.friction)*pipe.V_0**2 / (2*9.81*np.mean(pipe.diameter))
	pipe.attenuation = attenuation
	#print pipe.Name, pipe.friction.mean(),pipe.V_0,pipe.length,pipe.diameter,attenuation
	G.add_edge(pipe.node1,pipe.node2,tt = pipe.length / pipe.c, pipe = pipe,atten = attenuation,length = pipe.length)

	
ShortestPathLengths = nx.shortest_path_length(G,weight = 'tt')
ShortestPaths = nx.shortest_path(G,weight = 'tt')

Net.Graph = G

AttenuationArray = {}
for i in Net.nodes:
	for j in Net.nodes:
		
		TotalAtten = 1
		AttenuationArray[i,j] = TotalAtten
		Nodes = [ShortestPaths[i][j][0].Name]
		for k in range(1,len(ShortestPaths[i][j])):
			FrictionAtten = G.edges[ShortestPaths[i][j][k-1],ShortestPaths[i][j][k]]['atten']
			#Transmission =  ShortestPaths[i][j][k].TransCoeffs[G.edge[ShortestPaths[i][j][k-1]][ShortestPaths[i][j][k]]['pipe']]
			
			TotalAtten *= FrictionAtten
			AttenuationArray[i,j] *= FrictionAtten
			Nodes.append(ShortestPaths[i][j][k].Name)
		
		
		for k in range(1,len(ShortestPaths[i][j])-1):
			Transmission =  ShortestPaths[i][j][k].TransCoeffs[G.edges[ShortestPaths[i][j][k-1],ShortestPaths[i][j][k]]['pipe']]
			AttenuationArray[i,j] *= Transmission
			TotalAtten *= Transmission
			
		#print Nodes,TotalAtten
		j.locationMetric +=TotalAtten
	print i.Name
	
AttenArray2 = np.zeros((len(Net.nodes),len(Net.nodes))) 
iCount = 0
for i in Net.nodes:
	jCount = 0
	for j in Net.nodes:
		
		AttenArray2[iCount,jCount] = AttenuationArray[i,j]
		jCount+=1
	iCount +=1
	
	
pp.close('all')	
Net.best_Logger_Location(plot_Node_Names = False,No_to_highlight = 1)


	




SensorsRequired = 2
Best2 = np.zeros((len(Net.nodes),len(Net.nodes))) 
for i in range(len(Net.nodes)):
	for j in range(len(Net.nodes)):
		
		Best2[i,j] = sum(np.max(np.vstack((AttenArray2[:,i],AttenArray2[:,j])),axis=0))

ind = np.unravel_index(np.argmax(Best2, axis=None), Best2.shape)

node_highlight_list= []
for i in ind:
	node_highlight_list.append(Net.nodes[i])

Net.highlight_nodes(plot_Node_Names = False,nodes_highlighted = node_highlight_list)


#nx.draw_graphviz(G)
#Net.geom_Plot(plot_Node_Names = 1)	
#mng = pp.get_current_fig_manager()
#mng.window.wm_geometry("683x900+0+0")


#mng = pp.get_current_fig_manager()
#mng.window.wm_geometry("1366x900+0+0")
pp.show()
