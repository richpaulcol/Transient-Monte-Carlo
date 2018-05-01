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






Directory = 'Projects/WDSA_Models/'
FileName = '5Pipes.inp'
FileName = 'hanoi3.inp'
FileName = 'Net3LPSnoControl.inp'
Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()


#distribution = cp.J(cp.Normal(0.1,0.001))

#Iterations = 2
#samples = distribution.sample(Iterations)
#Output = np.zeros((4637,Iterations))

#for i in range(Iterations):
#	#Net.geom_Plot(plot_Node_Names = True)
#	Net.alter_epanet_friction(samples[i])
#	Net.run_epanet_file()
#	Net.read_results_from_epanet()
#	Wavespeed = 1159.2
#	Nodes = 20
#	Net.Constant_Wavespeed(Wavespeed)
#	dt = 50/(Wavespeed*Nodes)
#	Net.MOC_Initialisation(dt)


#	#Net.Assign_Emmiters('Emitters.csv')

##	Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
#	Net.Control_Input(Directory+'Demands.csv')
#	Net.MOC_Run(10)
#	Output[:,i] = Net.nodes[1].TranH
##	Output[0,i] = Net.nodes[1].H_0
#	
#Net.close_epanet_file()



#polynomial_expansion = cp.orth_ttr(3, distribution)
#foo_approx = cp.fit_regression(polynomial_expansion, samples[:10], Output[:,:10].T)

#expected = cp.E(foo_approx, distribution)
#deviation = cp.Std(foo_approx, distribution)

#x,y = normal_dist(expected[-1],deviation[-1]**2)
#x,y = normal_dist(expected,deviation**2)
