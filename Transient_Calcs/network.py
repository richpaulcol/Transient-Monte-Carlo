import numpy as np
import pylab as pp
from epanettools import epanet2 as epa
import time
#import networkx as nx
import scipy.sparse as ss
import pyhyd

""" 


Currently only Nodal Pressures are being saved



"""




def err(e):
    if(e>0):
        print e, epa.ENgeterror(e,25)
        #exit(5)


class Network(object):
	"""Represents a WDS Network."""
	def __init__(self,net_Name,nodes,pipes,valves,pumps,node_idx,pipe_idx,valve_idx,pump_idx):
		self.Name = net_Name
		self.nodes = nodes
		self.pipes = pipes
		self.valves = valves
		self.pumps = pumps
		
		#### Creating a list of link objects 
		self.links = self.pipes + self.valves + self.pumps
		
		self.node_idx = node_idx
		self.pipe_idx = pipe_idx
		self.valve_idx = valve_idx
		self.pump_idx = pump_idx
		
		#### Creating a dictionary of link objects id's
		self.link_idx = dict(self.pipe_idx.items() + self.valve_idx.items()+ self.pump_idx.items())
		
		
		#### Network Config Parameters
		self.Pressure_Dependent = 0
		self.epsilon = 0.85
		
		
		self.KE = [0]
		self.PE = [0]
		
		self.nodal_CPs = {}
		
	#####
	##	Function To Input EPANet Results as Initial Conditions
	
	def Import_EPANet_Results(self):	
		ret = epa.ENopen(self.filename,self.filename[:-3]+'rep',self.filename[:-3]+'out')
		print ret
		####	Opening the hydraulics results
		ret = epa.ENopenH()
		print ret
		ret = epa.ENinitH(0)
		print ret
		####	Running the Hydraulics Solver
		epa.ENrunH()
		
		####	Returning the head solutions (for the first time step, the only one accessible at the moment)
		##	Only need to do this for the node elements at the moment as reservoirs don't change
		ret,no_nodes=epa.ENgetcount(epa.EN_NODECOUNT)
		print 'Number of NODES in results file', no_nodes
		for index in range(1,no_nodes+1):
			ret,idx=epa.ENgetnodeid(index)
			ret,H0=epa.ENgetnodevalue(index, epa.EN_HEAD )
			try:
				#print Network.node_idx[idx].Name,idx
				if self.node_idx[idx].type == 'Node':
					self.node_idx[idx].H_0 = float(H0)
					self.node_idx[idx].TranH = [float(H0)]

			except:
				print 'Problem getting Head for Node:', idx
				continue
			
		####	Returning the flow solutions (for the first time step, the only one accessible at the moment)
		ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
		print 'Number of LINKS in results file',no_links
		for index in range(1,no_links+1):
			ret,idx=epa.ENgetlinkid(index)
			
			ret,Q0=epa.ENgetlinkvalue(index, epa.EN_FLOW )
			
			ret,V0=epa.ENgetlinkvalue(index, epa.EN_VELOCITY )

			ret,Headloss=epa.ENgetlinkvalue(index,epa.EN_HEADLOSS)
			ret,Length=epa.ENgetlinkvalue(index,epa.EN_LENGTH)
			ret,Diameter=epa.ENgetlinkvalue(index,epa.EN_DIAMETER)
			#print Headloss,Length,Diameter,V0
			#print 2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2)
			
			try:
			
				self.link_idx[idx].Q_0 = float(Q0)/1000. #Convert to m^3/s
				self.link_idx[idx].V_0 = float(V0)
				self.link_idx[idx].FF_0 = float(2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2))
			except:
				print 'Problem getting Flow or Velocity for link:', idx
				continue
			
	######
	##
	##	Initialise all links with a constant wave speed
	def Constant_Wavespeed(self,wavespeed):
		for i in self.links:
			i.c = wavespeed
			
	#######
	##	Initiailise all links with their respective wavespeeds
	def Initialise_Wavespeed(self):
		for i in self.links:
			i.Wavespeed()
	
	
	######
	##	Network MOC Initialisation
	##
	def MOC_Initialisation(self,dt):
		##  Start with a check of the dx values and the error in pipelengths and whether there is any steady state data
		self.dt = dt
		self.link_lengths = []
		self.link_dx = []
		self.link_length_error = []
		for i in self.links:
			try:
				i.Q_0
			except:
				print 'No steady-state data for link', i.Name
				break
			
			i.length = float(i.length)
			self.link_lengths.append(i.length)
			i.dx = i.c*dt
			self.link_dx.append(i.dx)
			i.NoTNodes = int(round(i.length / i.dx)+1)
			self.link_length_error.append(100.*abs((i.NoTNodes-1)*i.dx - i.length)/i.length)
			
			
			### Generating the transient constants
			#i.B = np.ones(i.NoTNodes)*i.c/(9.81*i.area)
			i.B = i.c/(9.81*i.area)
			
			### Generating the Transient Nodal positions
			i.TranNodesX = np.linspace(0,i.length,i.NoTNodes)
			### Initialising the Transient Node Heads and Flows
			
			"""The next two lines are removed to enable not saving the pipe pressure and flow"""
			#i.TranH = np.linspace(i.node1.H_0,i.node2.H_0,i.NoTNodes).reshape(1,i.NoTNodes)
			#i.TranQ = np.ones(i.TranNodesX.size).reshape(,i.NoTNodes)*i.Q_0
			
			i.TranH = np.linspace(i.node1.H_0,i.node2.H_0,i.NoTNodes).reshape(1,i.NoTNodes)
			i.TranH = np.vstack((i.TranH,i.TranH))
			i.TranQ = np.ones((2,i.TranNodesX.size))*i.Q_0
			
			
			### Setting up the pressure dependent demands for transient solutions
			if self.Pressure_Dependent == 1:
				for i in self.nodes:
					i.initialise_Pressure_Dependent_Demands()
					
		print 'Maximum length error =', max(self.link_length_error)
			
	#####
	##	Function to input the control file for the transient solver
	##	Currently only nodal flow control is available	
	##	Currently only works for numerical node names	
	def Control_Input(self,ControlFile):
		Control_Data = np.loadtxt(ControlFile,delimiter = ',')
		try:
			ControlNodes = Control_Data.shape[1]
		except:
			ControlNodes = 1
			Control_Data.resize(Control_Data.shape[0],1)

		for i in range(ControlNodes):
			self.node_idx[str(int(Control_Data[0,i]))].Transient_Demands = Control_Data[1:,i] / 1000. 
			self.node_idx[str(int(Control_Data[0,i]))].External_Flow = False ### Removes an emmiter from a control node, this is important otherwise you can't get flow rate control
		
			
	
		
	def Assign_Emmiters(self,EmitterFile):
		
		Emitter_Data = np.loadtxt(EmitterFile,delimiter = ',',dtype='string')

		
		Emitter_Nodes = Emitter_Data[0,1:]
		CdAs = Emitter_Data[1,1:]
		for i in range(0,len(CdAs)):
			try:
				CdAs[i] = float(CdAs[i])
			except:
				CdAs[i] = None
		
		#print Emitter_Nodes, CdAs
		for i in range(len(Emitter_Nodes)):
			#print Emitter_Nodes[i], CdAs[i], CdAs[i]=='None'
			self.node_idx[str(int(Emitter_Nodes[i]))].External_Flow = True
			if CdAs[i] == 'None':
				CdA = self.node_idx[str(int(Emitter_Nodes[i]))].demand / np.sqrt(self.node_idx[str(int(Emitter_Nodes[i]))].H_0)
				#print CdA
				self.node_idx[str(int(Emitter_Nodes[i]))].CdA = CdA 

			else:
				self.node_idx[str(int(Emitter_Nodes[i]))].CdA = float(CdAs[i])
	
	### Function to ensure a proper pressure based solution at all nodes that have a Demand
	### Sets an orifice at nodes with a demand with properties that match the SS demand.			
	def Assign_Emmiters_All(self):
		
		for i in range(len(self.nodes)):
			if self.nodes[i].type == 'Node':
				if self.nodes[i].demand != 0:				
					self.nodes[i].External_Flow = True
					self.nodes[i].CdA = self.nodes[i].demand / np.sqrt(self.nodes[i].H_0)
				
	def Assign_Random_Demands(self):
		for node in self.nodes:
			node.Random_Demand = True
	#####
	##	Function to run the transient solution
	##		
	def MOC_Run(self,maxT):
		self.Transient_Times = np.arange(0,maxT,self.dt)
		real_time = time.time()
		for t in self.Transient_Times[1:]:
			
			for i in self.links:
				#i.Name
				i.MOC_iter(self.epsilon,t)
			for i in self.nodes:
				#i.Name
				i.MOC_iter(self.epsilon,t,self.dt)	
			if t % 100*self.dt == 0:
				print 'Time:',t,'Elapsed Time',time.time() - real_time
				real_time = time.time()
		
		try:
			print 'Time:',t,'Elapsed Time',time.time() - real_time	
		except:	
			print 'No Time'
			
	
		
	######
	##	Plotting the network configuration
	def geom_Plot(self,plot_Node_Names = 0,Highlights = None):
		pp.figure()
		pp.title('Network Geometry')
		for i in self.nodes:
			if i.type == 'Node':
				symbol = 'o'
				size = 10
			elif i.type == 'Reservoir':
				symbol = 's'
				size = 20
			elif i.type == 'Tank':
				symbol = (5,2)
				size = 20	
			c = 'k'
			if i.Name in Highlights:
				size = 40	
				c = 'r'	
				
			pp.scatter([i.xPos],[i.yPos],marker = symbol,s = size,c=c)
			if plot_Node_Names != 0:
				pp.annotate(i.Name,(i.xPos,i.yPos))
			
		for i in self.pipes:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'k')
			
		for i in self.valves:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'r')
			
		for i in self.pumps:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'g')
		pp.axis('equal')
		pp.show()
		
		
	######
	##	PLotting the initial steady-state solution
	
	##	BROKEN
	
	
	def steady_Plot(self,Node_Range = None,Link_Range = None):
		Nodes_Data = []
		Links_Data = []
		for i in self.nodes:
			Nodes_Data.append([i.xPos,i.yPos,i.H_0])
		for i in self.links:
			Links_Data.append([[i.x1,i.x2],[i.y1,i.y2],i.Q_0])
		Nodes_Data = np.array(Nodes_Data).astype('float')
		Links_Data = np.array(Links_Data)
		pp.figure()
		
		if Node_Range == None:
			Node_max = Nodes_Data[:,2].max()
			
			Node_Colours = Nodes_Data[:,2] / Node_max
		else:
			Node_Colours = (Nodes_Data[:,2] - Node_Range[0]) / (Node_Range[1] - Node_Range[0])
		if Link_Range == None:
			Links_max = Links_Data[:,2].max()
			
			Link_Colours = Links_Data[:,2].astype('float') / Links_max
		else:
			Link_Colours = (Links_Data[:,2].astype('float') - Link_Range[0]) / (Link_Range[1] - Link_Range[0])	
		
		Link_Colours = Link_Colours.astype('float')
		pp.scatter(Nodes_Data[:,0],Nodes_Data[:,1],c = Node_Colours)
	#	for i in range(Links_Data.shape[1]):
	#		pp.plot(Links_Data[i,0],Links_Data[i,1],color = Link_Colours[i])
		
		pp.axis('equal')
		pp.show()
		pp.tight_layout()
		
		
	####
	##	Function to plot the transient solution at the network nodes.  If no additional argument is passed then all nodes are plotted.  The additional argument should consist of a list of node_ids
	def transient_Node_Plot(self, node_ids = 0):
		fig = pp.figure()
		ax = pp.subplot(111)
		pp.title('Nodal Head')
		if node_ids == 0:
			for i in self.node_idx:
				ax.plot(self.Transient_Times,self.node_idx[i].TranH,label = i)
		else:	
			for i in node_ids:
				ax.plot(self.Transient_Times,self.node_idx[i].TranH,label = i)
		ax.legend(loc='upper center',ncol=8)
		pp.xlabel('Time (s)')
		pp.ylabel('Head (m)')
		pp.show()
		pp.tight_layout()	
				
#	def best_Logger_Location(self,plot_Node_Names = 0,No_to_highlight = 0,No_Repeat_Length = 0):
#		locationMetrics = []
#		for i in self.nodes:
#			locationMetrics.append([i,i.locationMetric])
#		locationMetrics = np.array(locationMetrics)
#		#print locationMetrics
#		#locationMetrics.append(min(locationMetrics)/1.25)
#		pp.figure()
#		pp.title('Best Logger Locations')
#		for i in self.nodes:
#			if i.type == 'Node':
#				symbol = 'o'
#				size = 100
#			elif i.type == 'Reservoir':
#				symbol = 's'
#				size = 200
#			elif i.type == 'Tank':
#				symbol = (5,2)
#				size = 200	
#				
#			colour = (i.locationMetric - min(locationMetrics[:,1]))/( max(locationMetrics[:,1]) - min(locationMetrics[:,1]))
##			print colour
#			pp.scatter([i.xPos],[i.yPos],marker = symbol,s = size,c = 'r',alpha = colour)
#			if plot_Node_Names != 0:
#				pp.annotate(i.Name,(i.xPos,i.yPos))
#			
#		for i in self.pipes:
#			pp.plot([i.x1,i.x2],[i.y1,i.y2],'k')
#			
#		for i in self.valves:
#			pp.plot([i.x1,i.x2],[i.y1,i.y2],'r')
#			
#		for i in self.pumps:
#			pp.plot([i.x1,i.x2],[i.y1,i.y2],'g')
#		#print locationMetrics
#		Sorted = locationMetrics[locationMetrics.argsort(axis=0)[:,1],:]
#		#print Sorted
#		ShortestPathLengths = nx.shortest_path_length(self.Graph,weight = 'length')	
#		pp.scatter([Sorted[-1][0].xPos],[Sorted[-1][0].yPos],s=200, facecolors='none', edgecolors='k',linewidth = 4)
#		#print Sorted[-1][0].Name
#		for k in range(1,No_to_highlight):
#			#print Sorted[-k-1][0].Name
#			pp.scatter([Sorted[-k-1][0].xPos],[Sorted[-k-1][0].yPos],s=200, facecolors='none', edgecolors='k',linewidth = 4)
#			
#			
#			
#			#pp.annotate(str(k+1),(Sorted[-k-1][0].xPos,Sorted[-k-1][0].yPos))
#			
#		#if No_to_highlight >0:
#		#	self.Graph
#		
#			
#		pp.axis('equal')
#		pp.show()

	def set_fixed_friction(self,value):
		
		for pipe in self.pipes:
			pipe.FF_0 = value
			

	def open_epanet_file(self):
		ret = epa.ENopen(self.filename,self.filename[:-3]+'rep',self.filename[:-3]+'out')
	
	def run_epanet_file(self):
		####	Opening the hydraulics results
		ret = epa.ENopenH()
		print ret
		ret = epa.ENinitH(0)
		print ret
		####	Running the Hydraulics Solver
		epa.ENrunH()
		
	def read_results_from_epanet(self):
		####	Returning the head solutions (for the first time step, the only one accessible at the moment)
		ret,no_nodes=epa.ENgetcount(epa.EN_NODECOUNT)
		print 'Number of NODES in results file', no_nodes
		for index in range(1,no_nodes+1):
			ret,idx=epa.ENgetnodeid(index)
			ret,H0=epa.ENgetnodevalue(index, epa.EN_HEAD )
			ret,P0 = epa.ENgetnodevalue(index,epa.EN_PRESSURE)
			ret,Demand = epa.ENgetnodevalue(index, epa.EN_DEMAND)
			try:
				#print Network.node_idx[idx].Name,idx
				#if self.node_idx[idx].type == 'Node':
				self.node_idx[idx].H_0 = float(H0)
				self.node_idx[idx].P_0 = float(P0)
				self.node_idx[idx].TranH = [float(H0)]
				self.node_idx[idx].demand = float(Demand)/1000.

			except:
				print 'Problem getting Head for Node:', idx
				continue
			
		####	Returning the flow solutions (for the first time step, the only one accessible at the moment)
		ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
		print 'Number of LINKS in results file',no_links
		for index in range(1,no_links+1):
			ret,idx=epa.ENgetlinkid(index)
			
			ret,Q0=epa.ENgetlinkvalue(index, epa.EN_FLOW )
			
			ret,V0=epa.ENgetlinkvalue(index, epa.EN_VELOCITY )

			ret,Headloss=epa.ENgetlinkvalue(index,epa.EN_HEADLOSS)
			ret,Length=epa.ENgetlinkvalue(index,epa.EN_LENGTH)
			ret,Diameter=epa.ENgetlinkvalue(index,epa.EN_DIAMETER)
			ret,Roughness=epa.ENgetlinkvalue(index,epa.EN_ROUGHNESS)
			#print Headloss,Length,Diameter,V0
			#print 2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2)
			
			try:
			
				self.link_idx[idx].Q_0 = float(Q0)/1000. #Convert to m^3/s
				self.link_idx[idx].V_0 = float(V0)
				self.link_idx[idx].FF_0 = float(2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2))
				self.link_idx[idx].roughness = Roughness
				self.link_idx[idx].headloss = Headloss
			except:
				print 'Problem getting Flow or Velocity for link:', idx
				continue
				

	def alter_epanet_friction(self,value):
		ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
		for index in range(1,no_links+1):
			ret,epa.ENsetlinkvalue(index,epa.EN_ROUGHNESS,value[index-1])

	def alter_epanet_friction_all_same(self,value):
		ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
		for index in range(1,no_links+1):
			ret,epa.ENsetlinkvalue(index,epa.EN_ROUGHNESS,value)	

	def alter_epanet_demand(self,node,value):
		ret,index = epa.ENgetnodeindex(str(node))
		ret = epa.ENsetnodevalue(index,epa.EN_BASEDEMAND,value)
	
	def close_epanet_file(self):
		ret = epa.ENclose()
		
	def Initialise_Linear_Kalman(self,dt):
		self.dt = dt
		self.link_lengths = []
		self.link_dx = []
		self.link_length_error = []
		self.CPs = 0
		self.pipes_State_Index = np.array([])

		link_count = 0
		for i in self.links:
			try:
				i.Q_0
				
			except:
				print 'No steady-state data for link', i.Name
				break

			### Ensuring the friction calc is done
			#i.LambdaCalc()		#Usng the Friciton model as calculated from EPANet manual
			i.friction = i.FF_0 	#Using data back calculated from EPANet

			self.link_lengths.append(float(i.length))	##Length of Link
			i.dx = i.c*dt			
			self.link_dx.append(i.c*dt)					##Nodal distance for link
			i.NoCPs = int(round(i.length / i.dx)+1)			##Number of CPs for link
			self.link_length_error.append(100.*abs((i.NoCPs-1)*i.dx - i.length)/i.length)	## Max discretisation error
			i.B = i.c/(9.81*i.area)					## This is the basic calc constant 
			i.R = i.friction*i.dx*abs(i.Q_0) / (2*9.81*i.diameter*i.area**2)	## This is the basic friction cal which we can use if we assume that it is fixed throughout the calcualtions
			#i.R = i.R * i.dx
			i.TranCPsX = np.linspace(0,i.length,i.NoCPs)	##	Generating the CPs
			i.pipe_CP_index = (np.ones(i.TranCPsX.size)*link_count).astype('int')

			self.CPs += (i.TranCPsX).size	##working out the number of CPs
			self.pipes_State_Index = np.hstack((self.pipes_State_Index,i.pipe_CP_index))	##adding in the link to the state index
			
			i.CP_Node1 = self.pipes_State_Index.size - (i.TranCPsX).size
			i.CP_Node2 = self.pipes_State_Index.size -1
			i.link_number = link_count	##Adding in the link number to the link object
			link_count +=1
		
		self.X_Vector = (np.zeros((2*self.CPs+len(self.nodes)))).T #State Vector
		self.U_Vector = (np.zeros((2*self.CPs+len(self.nodes)))).T #Control Vector  (contains the demand values at each instance if not included in the state vector (which they could be and it would remain linear)  <- this should be step 2
		#print self.pipes_State_Index,self.pipes_State_Index.size
		

		####  	Generating the Matrices Required for the prediction steps
		self.A_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #Dynamics Evolution Matrix
		self.P_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #State Covariance Matrix
		self.Q_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #Uncertainty Covariance Matrix

		

		###	Looping again to generate the actual values in the X_Vector and A matrix
		for i in self.links:
			self.X_Vector[i.CP_Node1:i.CP_Node2+1] = np.linspace(i.node1.H_0,i.node2.H_0,i.NoCPs).reshape(1,i.NoCPs)	### The inital Head in the state vector
			self.X_Vector[self.CPs+i.CP_Node1:self.CPs+i.CP_Node2+1] = i.Q_0 #The initial Flow vectors
			
			#### The main sections along each pipe
			for CP in range(i.CP_Node1+1,i.CP_Node2):
				#print 'Central',CP,self.CPs,CP+self.CPs+1
				#Hp,Ha
				self.A_Matrix[CP,CP-1] = 0.5
				#Hp,Hb
				self.A_Matrix[CP,CP+1] = 0.5
				#Hp,Qa
				self.A_Matrix[CP,CP+self.CPs-1] = 0.5*i.B - 0.5*i.R
				#Hp,Qb
				self.A_Matrix[CP,CP+self.CPs+1] = -0.5*i.B + 0.5*i.R
				
				#Qp,Ha
				self.A_Matrix[CP+self.CPs,CP-1] = 1./(2.*i.B)
				#Qp,Hb
				self.A_Matrix[CP+self.CPs,CP+1] = -1./(2.*i.B)
				#Qp,Qa
				self.A_Matrix[CP+self.CPs,CP+self.CPs-1] = 0.5*(1-i.R/i.B)
				#Qp,Qb
				self.A_Matrix[CP+self.CPs,CP+self.CPs+1] = 0.5*(1-i.R/i.B)
		
		for i in self.nodes:
			
			###	Adding the CPs of the nodes to store the nodal head and flow values at the junctions
			try:
				i.nodal_CP = i.pipesIn[0].CP_Node2
			except:
				i.nodal_CP = i.pipesOut[0].CP_Node1
			self.nodal_CPs[i.idx] = i.nodal_CP
			#print i.type
			Bc = 0
			Qext = i.demand
			#print Qext,i.number
			for k in i.pipesOut:
				#print 'pOut',k.B
				Bc += 1./k.B
			for k in i.pipesIn:
				#print 'pIn',k.B
				Bc += 1./k.B
			Bc = 1./Bc
			#print Bc,'Pre'

			####
			##	Adding in the demand to the state vectors
			#print Qext
			self.X_Vector[2*self.CPs+i.number] = Qext
			self.A_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = 1
			#print self.X_Vector
			
			if i.type == 'Reservoir':
				#print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:
					#print 'Reservoir pipesIn',k.CP_Node2,k.CP_Node2+self.CPs,
					#Hp,Hp
					self.A_Matrix[k.CP_Node2,k.CP_Node2] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2] = -1./k.B
					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] = 1./k.B
					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] = 1. - k.R/k.B

				for k in i.pipesOut:
					#print 'Reservoir pipesOut',k.CP_Node1,k.CP_Node1+self.CPs,1. - k.R/k.B
					#Hp,Hp	
					self.A_Matrix[k.CP_Node1,k.CP_Node1] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1] = 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] = -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] = 1. - k.R/k.B

			if i.type == 'Node':
				#print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node2,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						self.A_Matrix[k.CP_Node2,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						self.A_Matrix[k.CP_Node2,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)
						#Qp,Ha
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node2,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						self.A_Matrix[k.CP_Node2,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						self.A_Matrix[k.CP_Node2,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)
						#Qp,Hb
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)

				for k in i.pipesOut:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node1,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						self.A_Matrix[k.CP_Node1,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						self.A_Matrix[k.CP_Node1,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)
						#Qp,Ha
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node1,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)
						#Qp,Hb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)

				for k in i.pipesIn:
					#Qp
					self.A_Matrix[k.CP_Node2+self.CPs,:] *= -1./k.B
#					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] += 1./k.B
##					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					self.A_Matrix[k.CP_Node2,2*self.CPs+i.number] = -Bc
					#Qa
					self.A_Matrix[k.CP_Node2+self.CPs,2*self.CPs+i.number] = Bc/k.B
					

				for k in i.pipesOut:
					#Qp
					self.A_Matrix[k.CP_Node1+self.CPs,:] *= 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] += -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					self.A_Matrix[k.CP_Node1,2*self.CPs+i.number] = -Bc
					#Qp
					self.A_Matrix[k.CP_Node1+self.CPs,2*self.CPs+i.number] = -Bc/k.B
		self.TransposeA = self.A_Matrix.T
		
		
	def regenerateA(self):
		Wavespeed = 1000.
		Diameter = 0.1
		Area = Diameter**2 *np.pi / 4
		Roughness = 1
		
		B = Wavespeed/(9.81*Area)					## This is the basic calc constant 
		
		Qs = self.X_Vector[self.CPs:2*self.CPs]
		Frictions = pyhyd.friction_factor(Diameter,Qs,Roughness/1000.,force_turb_flow =True)
		R = Frictions*self.dx*abs(Qs) / (2*9.81*Diameter*Area**2)
#		print R
		NewA_Matrix = np.zeros(self.A_Matrix.shape)
		NodePoints = self.nodal_CPs.values()
		NodePoints.sort()
		NodePoints[0] = -1
		for i in range(len(NodePoints)-1):
			for j in range(NodePoints[i]+2,NodePoints[i+1]):
				#Hp,Ha
				NewA_Matrix[j,j-1] = 0.5
				#Hp,Hb
				NewA_Matrix[j,j+1] = 0.5
				#Hp,Qa
				NewA_Matrix[j,j+self.CPs-1] = 0.5*B - 0.5*R[j-1]
				#Hp,Qb
				NewA_Matrix[j,j+self.CPs+1] = -0.5*B + 0.5*R[j+1]
				
				#Qp,Ha
				NewA_Matrix[j+self.CPs,j-1] = 1./(2.*B)
				#Qp,Hb
				NewA_Matrix[j+self.CPs,j+1] = -1./(2.*B)
				#Qp,Qa
				NewA_Matrix[j+self.CPs,j+self.CPs-1] = 0.5*(1-R[j-1]/B)
				#Qp,Qb
				NewA_Matrix[j+self.CPs,j+self.CPs+1] = 0.5*(1-R[j+1]/B)
		
		for i in self.nodes:
			## non pressure depedent
			#self.X_Vector[2*self.CPs+i.number] = Qext
			## pressure dependent
			#if i.demand != 0.:
			#	self.X_Vector[2*self.CPs+i.number] = np.sqrt(self.X_Vector[self.nodal_CPs[i.Name]]) * i.CdA
		
			NewA_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = 1   #Correct
			
			Bc= 0 
			for k in i.pipesOut:
				#print 'pOut',k.B
				Bc += 1./k.B
			for k in i.pipesIn:
				#print 'pIn',k.B
				Bc += 1./k.B
			Bc = 1./Bc
#			print Bc
			
			if i.type == 'Reservoir':
				#print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:
					#print 'Reservoir pipesIn',k.CP_Node2,k.CP_Node2+self.CPs,
					#Hp,Hp
					NewA_Matrix[k.CP_Node2,k.CP_Node2] = 1
					#Qp,Hp
					NewA_Matrix[k.CP_Node2+self.CPs,k.CP_Node2] = -1./k.B
					#Qp,Ha
					NewA_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] = 1./k.B
					#Qp,Qa
					NewA_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] = 1. - R[k.CP_Node2-1]/k.B

				for k in i.pipesOut:
					#print 'Reservoir pipesOut',k.CP_Node1,k.CP_Node1+self.CPs,1. - k.R/k.B
					#Hp,Hp	
					NewA_Matrix[k.CP_Node1,k.CP_Node1] = 1
					#Qp,Hp
					NewA_Matrix[k.CP_Node1+self.CPs,k.CP_Node1] = 1./k.B
					#Qp,Hb
					NewA_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] = -1./k.B
					#Qp,Qb
					NewA_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] = 1. - R[k.CP_Node1+1]/k.B

			if i.type == 'Node':
#				print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node2,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						NewA_Matrix[k.CP_Node2,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						NewA_Matrix[k.CP_Node2,j.CP_Node2+self.CPs-1] = Bc*(1-R[j.CP_Node2-1]/j.B)
#						#Qp,Ha
						NewA_Matrix[k.CP_Node2+self.CPs,j.CP_Node2-1] = Bc/j.B
#						Qp,Qa
						NewA_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-R[j.CP_Node2-1]/j.B)
#						NewA_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node2,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						NewA_Matrix[k.CP_Node2,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						NewA_Matrix[k.CP_Node2,j.CP_Node1+self.CPs+1] = -Bc*(1-R[j.CP_Node1+1]/j.B)
#						#Qp,Hb
						NewA_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						NewA_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-R[j.CP_Node1+1]/j.B)

				for k in i.pipesOut:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node1,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						NewA_Matrix[k.CP_Node1,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						NewA_Matrix[k.CP_Node1,j.CP_Node2+self.CPs-1] = Bc*(1-R[j.CP_Node2-1]/j.B)
						#Qp,Ha
						NewA_Matrix[k.CP_Node1+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						NewA_Matrix[k.CP_Node1+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-R[j.CP_Node2-1]/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node1,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						NewA_Matrix[k.CP_Node1,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						NewA_Matrix[k.CP_Node1,j.CP_Node1+self.CPs+1] = -Bc*(1-R[j.CP_Node1+1]/j.B)
						#Qp,Hb
						NewA_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						NewA_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-R[j.CP_Node1+1]/j.B)

				for k in i.pipesIn:
					#print Bc,k.B,Bc * Qext / k.B
					#print 'Dave'
					#print 'Final Node pipesIn',k.CP_Node2,k.CP_Node2+self.CPs,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
					#self.U_Vector[k.CP_Node2] = -Bc * Qext #Adding the control to the head node				
					#self.U_Vector[k.CP_Node2+self.CPs] = Bc * Qext /k.B#Adding the control 

					
					#Qp
					NewA_Matrix[k.CP_Node2+self.CPs,:] *= -1./k.B
#					#Qp,Ha
					NewA_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] += 1./k.B
##					#Qp,Qa
					NewA_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					NewA_Matrix[k.CP_Node2,2*self.CPs+i.number] = -Bc
					#Qa
					NewA_Matrix[k.CP_Node2+self.CPs,2*self.CPs+i.number] = Bc/k.B
#					

				for k in i.pipesOut:
					#print Bc,k.B,-Bc * Qext / k.B
					#print 'Final Node pipesOut',k.CP_Node1,k.CP_Node1+self.CPs,j.CP_Node1+1,j.CP_Node1+self.CPs+1
					#self.U_Vector[k.CP_Node1] = -Bc* Qext #Adding the control to the head node				
					#self.U_Vector[k.CP_Node1+self.CPs] = -Bc * Qext / k.B #Adding the control 
					
					
					#Qp
					NewA_Matrix[k.CP_Node1+self.CPs,:] *= 1./k.B
					#Qp,Hb
					NewA_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] += -1./k.B
					#Qp,Qb
					NewA_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					NewA_Matrix[k.CP_Node1,2*self.CPs+i.number] = -Bc
					#Qp
					NewA_Matrix[k.CP_Node1+self.CPs,2*self.CPs+i.number] = -Bc/k.B
#			
		return NewA_Matrix
		
		
	def UnscentedInitialise(self,dt,InputVariables):
		self.dt = dt
		self.link_lengths = []
		self.link_dx = []
		self.link_length_error = []
		self.CPs = 0
		self.pipes_State_Index = np.array([])
		self.InputVariables = InputVariables

		link_count = 0
		for i in self.links:
			try:
				i.Q_0
				
			except:
				print 'No steady-state data for link', i.Name
				break

			### Ensuring the friction calc is done
			#i.LambdaCalc()		#Usng the Friciton model as calculated from EPANet manual
			i.friction = i.FF_0 	#Using data back calculated from EPANet

			self.link_lengths.append(float(i.length))	##Length of Link
			i.dx = i.c*dt			
			self.link_dx.append(i.c*dt)					##Nodal distance for link
			i.NoCPs = int(round(i.length / i.dx)+1)			##Number of CPs for link
			
			self.link_length_error.append(100.*abs((i.NoCPs-1)*i.dx - i.length)/i.length)	## Max discretisation error

			i.B = i.c/(9.81*i.area)					## This is the basic calc constant 
			i.R = i.friction*i.dx*abs(i.Q_0) / (2*9.81*i.diameter*i.area**2)	## This is the basic friction cal which we can use if we assume that it is fixed throughout the calcualtions
			#i.R = i.R * i.dx
			i.TranCPsX = np.linspace(0,i.length,i.NoCPs)	##	Generating the CPs
			i.pipe_CP_index = (np.ones(i.TranCPsX.size)*link_count).astype('int')

			self.CPs += (i.TranCPsX).size	##working out the number of CPs
			self.pipes_State_Index = np.hstack((self.pipes_State_Index,i.pipe_CP_index))	##adding in the link to the state index
			
			i.CP_Node1 = self.pipes_State_Index.size - (i.TranCPsX).size
			i.CP_Node2 = self.pipes_State_Index.size -1
			i.link_number = link_count	##Adding in the link number to the link object
			link_count +=1
		
		self.X_Vector = (np.zeros((2*self.CPs+InputVariables))).T 
		
	def MOCtoStateVector(self,State):
		for i in range(len(self.pipes)):
			State[:self.CPs][self.pipes_State_Index==i] = self.pipes[i].TranH[1,:]
			State[self.CPs:2*self.CPs][self.pipes_State_Index==i] = self.pipes[i].TranQ[1,:]
		return State
	
	
	
	def StateVectortoMOC(self,State):
		### Need to put the info into the TranH and TranQ [0,:] vector
		for i in range(len(self.pipes)):
			self.pipes[i].TranH[0,:] = State[:self.CPs][self.pipes_State_Index==i]
			self.pipes[i].TranQ[0,:] = State[self.CPs:2*self.CPs][self.pipes_State_Index==i]
			self.pipes[i].roughness = np.abs(State[2*self.CPs+i])
		
		
		self.nodes[-1].H_0  = State[0]
		
		self.nodes[1].demand = State[-1]
		self.nodes[4].demand = State[-2]
		#self.pipes[0].TranH[0,0] = State[-3]
		
		

	
	def UpdateState(self,State,noise = 0):
		#print noise
		self.StateVectortoMOC(State+noise)
		
		for i in self.links:
			i.MOC_iter(self.epsilon,self.time)
		for i in self.nodes:		
			i.MOC_iter(self.epsilon,self.time,self.dt)		
				
		State = self.MOCtoStateVector(State)
		return State
