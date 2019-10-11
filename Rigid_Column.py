#import numpy as np
#import pylab as pp

#import wntr as wntr
#import networkx as nx

""" This file at the moment implements the explicit incidence method for rwc calculations from Shimada 1989.

It currently doesn't have pressure dependent demands, or a true non-linear friction model, nor does it allow for tanks.

I am also using the notation from Nault (Thesis) 2016.

TODO I need to explore the stability of the code and when the explicit formulation will break down,
but for now it is easier to code and understand."""

from Transient_Calcs import *
from numpy.linalg import multi_dot

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

	Net.A = np.zeros((NodesNum, LinksNum))
	Net.A0 = np.zeros((len(Net.Type0), LinksNum))
	Net.A1 = np.zeros((len(Net.Type1), LinksNum))
	for i in Net.links:
		LinkNumber = i.number
		LinkStart = i.node1.number
		LinkEnd = i.node2.number
		Net.A[LinkStart, LinkNumber] = 1
		Net.A[LinkEnd, LinkNumber] = -1

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
	Net.Q = np.ones((LinksNum))
	Net.q0 = np.zeros((len(Net.Type0)))
	Net.q1 = np.zeros((len(Net.Type1)))
	Net.dq1 = np.zeros((len(Net.Type1)))
	Net.H0 = np.zeros((len(Net.Type0)))
	Net.H1 = np.zeros((len(Net.Type1)))
	Net.L = np.zeros((LinksNum,LinksNum))
	Net.InvL = np.zeros((LinksNum,LinksNum))
	Net.C1 = np.zeros((len(Net.Type1),len(Net.Type1)))
	Net.InvC1 = np.zeros((len(Net.Type1),len(Net.Type1)))
	Net.F = np.zeros((LinksNum))

def calc_F_submatrix(Net):
	LinksNum = len(Net.links)
	for i in range(LinksNum):
		if Net.links[i].linktype == 'Pipe':
			Net.links[i].Re = (abs(Net.links[i].Q_0)/(Net.links[i].area) * Net.links[i].diameter / (1.004*10.**-6))
			Net.links[i].Friction()
			#print Net.links[i].friction,Net.links[i].FF_0


			Net.F[i] = 8*Net.links[i].FF_0 * Net.links[i].length / (g*np.pi*Net.links[i].diameter**5) * Net.Q[i]*abs(Net.Q[i])

		elif Net.links[i].linktype == 'Valve':
			E = 1.0
			Net.F[i] = 1./(Net.links[i].setting**2 * E**2) * Net.Q[i]*abs(Net.Q[i])

def create_Jacobian_submatrix(Net):
	LinksNum = len(Net.links)
	Net.D = np.zeros((LinksNum,LinksNum))
	if Net.Friction_Units == 'H-W':
		n = 1.852

		for i in range(LinksNum):
			#if Net.links[i].type == 'Pipe'
			K = 10.67*Net.links[i].length / (Net.links[i].diameter**4.87 * Net.links[i].roughness)
			Net.D[i,i] = n * K * abs(Net.Q[i])**(n-1)


	if Net.Friction_Units == 'D-W':
		n=2
		for i in range(LinksNum):
			K = 8*Net.links[i].length*Net.links[i].friction / (9.81 * np.pi * Net.links[i].diameter**5)
			Net.D[i,i] = n * K * abs(Net.Q[i])**(n-1)

def create_M(Net):
	LinksNum = len(Net.links)
	Net.M = np.zeros((LinksNum, LinksNum))
	Net.M_1 = np.zeros((LinksNum, LinksNum))
	for i in range(LinksNum):
		Net.M[i,i] = Net.links[i].length / (9.81 * Net.links[i].area) + Net.D[i,i]
		Net.M_1[i,i] = 1./Net.M[i,i]

def rwc_iteration(Net,Q,H,q,dt):
	#Todo Produce a vectorised function here for a RWC iteration.
	print 'Dave'


Directory = 'Projects/Water_Seminar/'
FileName = 'Net3b_no_Tank.LPS.inp'

Directory = 'Projects/SimpleBranched/'
FileName = 'SimpleBranchedTrue.inp'
FileName = '5_Pipes.inp'

Directory = 'Projects/Transient_Fingerprinting/'
FileName = 'transient vlbv4x cut down.inp'


Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()
create_Incidence_Matrix(Net)
create_Solution_Vectors(Net)
# create_Jacobian_submatrix(Net)
Net.run_epanet_file()
Net.read_results_from_epanet()

for i in range(len(Net.links)):
	Net.Q[i] = Net.links[i].Q_0
	Net.L[i,i] = Net.links[i].length / (9.81 * Net.links[i].area)
	Net.InvL[i,i] = 1./Net.L[i,i]

Type0count = 0
Type1count = 0
for i in range(len(Net.nodes)):
	if Net.nodes[i].type == 'Reservoir':
		Net.H0[Type0count] = Net.nodes[i].H_0
		Net.q0[Type0count] = Net.nodes[i].demand
		Type0count +=1
	else:
		Net.H1[Type1count] = Net.nodes[i].H_0
		Net.q1[Type1count] = Net.nodes[i].demand
		Net.C1[Type1count,Type1count] = 0.0001
		Net.InvC1[Type1count,Type1count] = 1./Net.C1[Type1count,Type1count]
		Type1count +=1

calc_F_submatrix(Net)

### Explicit Integration

Qs = []#Net.Q]
Hs = []#Net.H1]
qs = []



#  My method, (Doesn't seem to work)
# for iters in range(100000):
# 	calc_F_submatrix(Net)
#
# 	dQ = multi_dot((-Net.InvL,Net.F))+multi_dot((Net.InvL,Net.A0.T,Net.H0)) + multi_dot((Net.InvL,Net.A1.T,Net.H1))
# 	dH = multi_dot((-Net.InvC1,Net.A1,Net.Q)) - multi_dot((-Net.InvC1,Net.q1))
#
# 	Net.Q = Net.Q + dt*dQ
# 	Net.H1 = Net.H1 + dt*dH
# 	Qs.append(Net.Q)
# 	Hs.append(Net.H1)
#
# Qs = np.array(Qs)
# Hs = np.array(Hs)
#
# f,axs = pp.subplots(2,1)
# axs[0].plot(Qs)
# axs[1].plot(Hs)
# pp.show()

Net.I = np.identity(len(Net.links))
Net.R = multi_dot((Net.A1,Net.InvL,Net.A1.T))
Net.InvR = np.linalg.inv(Net.R)
Net.S = multi_dot((Net.A1.T,Net.InvR,Net.A1,Net.InvL))
Net.W = multi_dot((Net.InvL,(Net.I - Net.S)))

dt = 1./200.
maxt = 1

times = np.arange(0,maxt,dt)
start_Time = 5
start_Iter = int(start_Time / dt)
end_Time = 55
end_Iter = int(end_Time / dt)

ControlNode = 1
Control = np.ones(times.shape)*Net.q1[ControlNode]
Control[start_Iter:] = (0.6+0.4*np.cos(1*times[:-start_Iter]))*Net.q1[ControlNode]
#Control[Start_Iter:End_Iter] =
#
# Control2 = np.zeros(times.shape)
# Control2[start_Iter:] = 0.001



for time in times:


	Net.dq1[ControlNode] = (pp.gradient(Control)/dt)[int(time/dt)]
	Net.q1[ControlNode] = Control[int(time/dt)]



	#Net.dq1[9] = (pp.gradient(Control2)/dt)[int(time/dt)]
	#Net.q1[9] = Control2[int(time / dt)]

	if time >= 0.5:
		Net.valves[2].setting = 1e-6
	# 	Net.dq1[8] = -0.005 / dt
	# else:
	# 	Net.dq1[8] = 0

	#if time > 2:
#		Net.H0[0] = 300

	#Net.q1 = Net.q1 + Net.dq1

	calc_F_submatrix(Net)
	T = multi_dot((Net.InvL, Net.A1.T, Net.InvR, Net.dq1))
	dQ = multi_dot((-Net.W, Net.F)) + multi_dot((Net.W, Net.A0.T, Net.H0))-T
	dQ = -T
	Net.OldQ = Net.Q
	Net.Oldq1 = Net.q1
	Net.Q = Net.Q + dt * dQ
	Net.H1 =multi_dot((Net.InvR,multi_dot((Net.A1,Net.InvL,(Net.F - multi_dot((Net.A0.T,Net.H0))))) - Net.dq1))
	#Net.newq1 =
	#Net.dq1 = (Net.newq1-Net.q1)/dt
	Net.q1 = multi_dot((Net.A1,Net.Q))
	Net.q0 = multi_dot((Net.A0,Net.Q))

	Net.q = np.hstack((Net.q0,Net.q1))
	Net.H = np.hstack((Net.H0,Net.H1))

	Qs.append(Net.Q)
	Hs.append(Net.H)
	qs.append(Net.q)
#
Qs = np.array(Qs)
Hs = np.array(Hs)
qs = np.array(qs)

f,axs = pp.subplots(3,1)
axs[0].plot(times,Qs)
axs[1].plot(times,Hs)
axs[2].plot(times,qs)
pp.show()


#Todo: I need to recreate this file run in the standard MOC and see how closely the two results are.
#Todo: Implement this RWC in a unscented Kalman filter
