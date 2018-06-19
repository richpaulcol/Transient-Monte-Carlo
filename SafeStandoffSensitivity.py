import pylab as pp
import numpy as np
from Transient_Calcs import *


Directory = 'Projects/SafeStandOff/'
FileName = 'Sensitivity.inp'


Net = Import_EPANet_Geom(Directory+FileName)
Net.open_epanet_file()

Net.run_epanet_file()
Net.read_results_from_epanet()

PipeDiameterRange = np.arange(0.1,1.6,0.1)
OrificeDiameterRange = np.arange(0.,1.0,0.1)
ReservoirHeadRange = np.arange(0,100,20)
OrificeCoefficientRange = np.arange(0.5,0.9,0.05)
PipeLengthRange = np.arange(100,1000,100)

Inputs = np.meshgrid(PipeDiameterRange,OrificeDiameterRange,ReservoirHeadRange,OrificeCoefficientRange,PipeLengthRange,indexing ='ij')
WSSC = np.zeros(Inputs[0].shape)
DUTCH1 = np.zeros(Inputs[0].shape)
DUTCH2 = np.zeros(Inputs[0].shape)

node = Net.nodes[0]
ret,node_index = epa.ENgetnodeindex(str(node.Name))
pipe = Net.pipes[0]
ret,pipe_index = epa.ENgetlinkindex(str(pipe.Name))
res = Net.nodes[1]
ret,res_index = epa.ENgetnodeindex(str(res.Name))

for i in range(PipeDiameterRange.size):
	for j in range(OrificeDiameterRange.size):
		for k in range(ReservoirHeadRange.size):
			for l in range(OrificeCoefficientRange.size):
				for m in range(PipeLengthRange.size):
					AO = Inputs[1][i,j,k,l,m]**2*np.pi/4	
					
					ret = epa.ENsetnodevalue(node_index,epa.EN_EMITTER,1000*Inputs[3][i,j,k,l,m]*AO*np.sqrt(2*9.81))
					
					ret = epa.ENsetnodevalue(res_index,epa.EN_ELEVATION,Inputs[2][i,j,k,l,m])
					
					ret = epa.ENsetlinkvalue(pipe_index,epa.EN_DIAMETER,Inputs[0][i,j,k,l,m])
					ret = epa.ENsetlinkvalue(pipe_index,epa.EN_LENGTH,Inputs[4][i,j,k,l,m])

					Net.run_epanet_file()
					Net.read_results_from_epanet()

					ret,nodeQ = epa.ENgetnodevalue(node_index,epa.EN_DEMAND)
					ret,nodeH = epa.ENgetnodevalue(node_index,epa.EN_HEAD)
					print(1000*Inputs[3][i,j,k,l,m]*AO*np.sqrt(2*9.81),nodeQ,nodeH)

					WSSC[i,j,k,l,m] = (1/1.5)*(Inputs[3][i,j,k,l,m]**2*Inputs[2][i,j,k,l,m])*(1-np.tan(np.radians(30)))
					DUTCH1[i,j,k,l,m] = 8.0 * (Inputs[2][i,j,k,l,m]**3 * Inputs[0][i,j,k,l,m]**5)**(1./8.)
					DUTCH1[i,j,k,l,m] = 1.8 *(nodeQ * nodeH/(9.81**0.5 * Inputs[3][i,j,k,l,m] * Inputs[1][i,j,k,l,m]**(7./2.)))**0.243 * Inputs[1][i,j,k,l,m]
epa.ENsaveinpfile(Directory+'Final.inp')
