import pylab as pp
import numpy as np
import chaospy as cp
import os
from Transient_Calcs import *
import pykalman as pk
#import seaborn as sns

from numpy import linalg as la

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3
    
def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

def rescaleCOV(COV,scale):
	from scipy import interpolate
	mymin,mymax = 0,COV.shape[0]-1
	X = np.linspace(mymin,mymax,COV.shape[0])
	Y = np.linspace(mymin,mymax,COV.shape[0])
	x,y = np.meshgrid(X,Y)
	f = interpolate.interp2d(x,y,COV,kind='linear')
	Xnew = np.linspace(mymin,mymax,(COV.shape[0]-1)*scale)                      
	Ynew = np.linspace(mymin,mymax,(COV.shape[0]-1)*scale)
	return (f(Xnew,Ynew)+f(Xnew,Ynew).T)/2.

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

Directory = 'Projects/Kalman_Pipe_Paper/'
FileName = 'SinglePipe.inp'
maxTime = 30
dt =0.01
Wavespeed = 1000.
Transient_Times = np.arange(0,maxTime,dt)

#####
# 	The input distributions
#	Upstream Head =  N(20,0.1)  [m]
#	Downstream Flow = N(0.001,0.0001)   [m3/s]
#	Demand at node3 = Exp(0.0005)
#	Roughness = N(1,0.4)

#Upstream = cp.Normal(20,0.1)
#Downstream = cp.Normal(0.001,0.0001)
#Demand = cp.Exponential(0.0005)
#Roughness = cp.J(cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4))

#Distributions = cp.J(Upstream,Downstream,Demand,Roughness)
#samples = Distributions.sample(100000,'S')
#np.save(Directory+'Transient_MC_samples.npy',samples)


#

#####		Linear Transient
#
#	Calculating the mean inputs

Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()

#####	Setting up the initial condition
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

######	Running the initial steady state
Net.run_epanet_file()				## This function runs the SS EPAnet model	
Net.read_results_from_epanet()			## Reads the EPAnet data
Net.Constant_Wavespeed(Wavespeed)		## Sets all pipes with the same wavespeed
Net.Initialise_Linear_Kalman(dt)		## Descritises the model and sets up all the matrices
Net.Assign_Emmiters_All()			## Sets all SS outputs to be emmiters rather than constant flow
Net.dx = Wavespeed*dt				## Calculates a 
#####	Initialising the Covariance Matrices
#COV = np.load(Directory+'InitialCOV.npy')
#Var = np.diag(COV)

Net.P_Matrix = Net.P_Matrix.todense()
Net.Q_Matrix = Net.Q_Matrix.todense()
##for i in range(Net.pipes_State_Index.size):
##	for j in range(Net.pipes_State_Index.size):
##		k = Net.pipes_State_Index.astype('int')[i]
##		l = Net.pipes_State_Index.astype('int')[j]
##		Net.P_Matrix[Net.CPs+i,Net.CPs+j] = COV[k,l]
##		Net.P_Matrix[i,j] = COV[k+len(Net.pipes),l+len(Net.pipes)]
#		
#		
##Net.P_Matrix[:Net.CPs,:Net.CPs] = (np.ones((50,50))*0.15)
##Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs] = (np.ones((50,50))*0.01)
##Net.P_Matrix[Net.CPs:Net.CPs+20,Net.CPs:Net.CPs+20] = (np.ones((20,20))*0.26)
#Net.P_Matrix[0,0] = 0.1**2  				#Variance in the upstream Head BC#  (i.e. std^2)
#Net.P_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.0001**2	#Variance in the downstream flow#
#Net.P_Matrix[2*Net.CPs+1,2*Net.CPs+1] = 0.0005**2	#Variance in the node 3 demand

#Net.P_Matrix[:Net.CPs,:Net.CPs] = rescaleCOV(COV[5:,5:],10)

#Net.P_Matrix = np.load(Directory+'InitP.npy')
##Net.Q_Matrix[0,0] = 0.1**2*Net.dt**2  		#Variance in the upstream Head BC#
##Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = Net.dt**2*0.0001**2	#Variance in the downstream flow#
##Net.Q_Matrix[2*Net.CPs+1,2*Net.CPs+1] = Net.dt**2*0.0005**2	#Variance in the node 3 demand


#sigma_lam = 0.00425


##Net.Q_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs] = np.diag(8 * Net.dx* Net.X_Vector[Net.CPs:2*Net.CPs]**2 / (9.81*0.1**5 * np.pi**2)) * sigma_lam**2


##Net.Q_Matrix[Net.CPs,Net.CPs] = (4*Net.dx)/(9.81*0.1**5*np.pi)*(-Net.X_Vector[Net.CPs+1]**2 + Net.X_Vector[Net.CPs]**2)*sigma_lam**2
##Net.Q_Matrix[2*Net.CPs-1,2*Net.CPs-1] = (4*Net.dx)/(9.81*0.1**5*np.pi)*(-Net.X_Vector[2*Net.CPs-1]**2 + Net.X_Vector[2*Net.CPs-2]**2)*sigma_lam**2



####### Adding the uncertainty to the head on the pipes due to friction 
##Net.Q_Matrix[1:Net.CPs-1,1:Net.CPs-1] =  (4*Net.dx)/(9.81*0.1**5*np.pi)  *  (-Net.X_Vector[Net.CPs+2:2*Net.CPs]**2 + Net.X_Vector[Net.CPs:2*Net.CPs-2]**2)*sigma_lam**2* Net.dt/2.


####### Adding the uncertainty to the flow on the pipes due to friction 
##Net.Q_Matrix[Net.CPs+1:2*Net.CPs-1,Net.CPs+1:2*Net.CPs-1] = (Net.dx/(0.1**3 * np.pi*Wavespeed))  *  (+Net.X_Vector[Net.CPs+2:2*Net.CPs]**2 + Net.X_Vector[Net.CPs:2*Net.CPs-2]**2)  *sigma_lam**2 * Net.dt/2.

##fsfsdfsfs


#kf = pk.KalmanFilter()

#kf.transition_matrices = Net.A_Matrix.todense()

#kf.transition_covariance = Net.Q_Matrix

#kf.initial_state_mean = Net.X_Vector
#kf.initial_state_covariance = Net.P_Matrix
#Simon = kf.transition_matrices
#Dave = Net.regenerateA()
#StateOutput = np.zeros((Net.X_Vector.size,Transient_Times.size))
#VarianceOutput = np.zeros((Net.X_Vector.size,Net.X_Vector.size,Transient_Times.size))
#for i in range(0,Transient_Times.size):
#	if i == int(2/dt):
#		Net.X_Vector[2*Net.CPs+Net.nodes[-2].number]+=0.0005
#		#Net.P_Matrix[2*Net.CPs+4,2*Net.CPs+4] += 0.00001**2
#		
#	Net.X_Vector,Net.P_Matrix = kf.filter_update(Net.X_Vector,Net.P_Matrix)
#	StateOutput[:,i] = Net.X_Vector
#	VarianceOutput[:,:,i] = Net.P_Matrix 


#	kf.transition_matrices = Net.regenerateA()
##	Net.Q_Matrix[1:Net.CPs-1,1:Net.CPs-1] =  (4*Net.dx)/(9.81*0.1**5*np.pi)  *  (-Net.X_Vector[Net.CPs+2:2*Net.CPs]**2 + Net.X_Vector[Net.CPs:2*Net.CPs-2]**2)*sigma_lam**2* Net.dt/2.
##	Net.Q_Matrix[Net.CPs+1:2*Net.CPs-1,Net.CPs+1:2*Net.CPs-1] = (Net.dx/(0.1**3 * np.pi*Wavespeed))  *  (+Net.X_Vector[Net.CPs+2:2*Net.CPs]**2 + Net.X_Vector[Net.CPs:2*Net.CPs-2]**2)  *sigma_lam**2 * Net.dt/2.

#np.save(Directory+'KalmanStateOutput.npy',StateOutput)
#np.save(Directory+'KalmanVarianceOutput.npy',VarianceOutput)

#import Pressure_Surges_Monte_Plotter
##f,axs = pp.subplots(nrows = 2,ncols = 1)
#axs[0].plot(Transient_Times,StateOutput[Net.nodal_CPs['3'],:])
#axs[1].plot(Transient_Times,VarianceOutput[Net.nodal_CPs['3'],Net.nodal_CPs['3'],:]*2)

#f,axs = pp.subplots(nrows = 4,ncols = 2)
#axs[0,0].imshow(Simon[:50,:50])
#axs[0,1].imshow(kf.transition_matrices[:50,:50])

#axs[1,0].imshow(Simon[50:,50:])
#axs[1,1].imshow(kf.transition_matrices[50:,50:])

#axs[2,0].imshow(Simon[:50,50:])
#axs[2,1].imshow(kf.transition_matrices[:50,50:])

#axs[3,0].imshow(Simon[50:,:50])
#axs[3,1].imshow(kf.transition_matrices[50:,:50])

#pp.show()



