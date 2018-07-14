import numpy as np
import pylab as pp
import seaborn as sns

def normal_dist(mean,variance):
	stddev = np.sqrt(variance)
	x = np.linspace(mean-3*stddev,mean+3*stddev,100)
	y = 1./(np.sqrt(2 * np.pi * variance)) * np.exp(- (x-mean)**2 / (2*variance))
	return x,y


Directory = 'Projects/Pressure_Surges_Models/'
FileName = '5Pipes.inp'


Outputs = np.load(Directory+'InitialConditions.npy')
Outputs[:,:5] = Outputs[:,:5]/1000.

### Plot of the correlation and covariance matrices

COV = np.cov(Outputs.T)
np.save(Directory+'InitialCOV.npy',COV)
COR = np.corrcoef(Outputs.T)

AxisLabels = ['Q1','Q2','Q3','Q4','Q5','H1','H2','H3','H4','H5','H6']
cmap = pp.get_cmap('jet', 30)
f, ax = pp.subplots(figsize=(9, 6))                        
sns.heatmap(COV, fmt="0.3f",annot=True, linewidths=.5, ax=ax,cmap=cmap,xticklabels = AxisLabels,yticklabels=AxisLabels)    
#sns.heatmap(COR, fmt="0.3f",annot=True, linewidths=.5, ax=ax,cmap=cmap,xticklabels = AxisLabels,yticklabels=AxisLabels)              
pp.tight_layout()
pp.show()
pp.savefig(Directory+'5_Pipes_Initial_Covariance_Matrix.png')

### Plot of the convergence of the mean and variance
#means = []
#variances = []
#NoSamples = np.arange(2,6.1,0.1)
##NoSamples = [10,30,100,300,1000,3000,10000,30000,100000,300000,1000000]
#for i in 10**NoSamples:
#	means.append(np.mean(Outputs[:int(i),-1]))
#	variances.append(np.std(Outputs[:int(i),-1])**2)
#	

#pp.figure(figsize=(6, 4))
#pp.semilogx(10**NoSamples,means)
#pp.xlabel('Monte Carlo Samples')
#pp.ylabel('Mean of Head Node 6')
#pp.tight_layout()
#pp.savefig(Directory+'5_Pipes_Initial_Convergence_Mean.png')
#pp.show()


