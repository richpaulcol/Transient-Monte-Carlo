import pylab as pp
import numpy as np
import chaospy as cp
import os


Directory = 'Projects/Pressure_Surges_Models/'
FileName = '5Pipes.inp'
maxTime = 30
dt =0.01
Transient_Times = np.arange(0,maxTime,dt)

distribution = cp.J(cp.Normal(1,0.4))#,cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4),cp.Normal(1,0.4))
demand_distribution = cp.Exponential(1.)

Iterations = 1000

#samples = np.load(Directory + '5_pipes_varying_demand_samples.npy')
#output = np.load(Directory + '5_pipes_varying_demand.npy')

samples = np.load(Directory + '5_pipes_varying_friction_all_same_samples.npy')
output = np.load(Directory + '5_pipes_varying_friction_all_same.npy')

####
#	Plotting input distributions
#pp.figure()
#pp.hist(samples[0,:],density=True,color='k',bins=40,alpha=0.50,label='Samples')
#d = np.linspace(min(samples[0,:]),max(samples[0,:]),1000)
#pp.plot(d,cp.Normal(1,0.4).pdf(d),'k',label='Roughness PDF')
#pp.xlabel('Roughness (mm)')
#pp.ylabel('Frequency')
#pp.show()

Node = 4
ConditionPoints = 100
Order = 3

polynomial_expansion = cp.orth_ttr(Order, distribution)
foo_approx = cp.fit_regression(polynomial_expansion, samples[:ConditionPoints], output[:ConditionPoints,Node,:])

coefs_kernal = cp.descriptives.misc.QoI_Dist(foo_approx,distribution)


#pp.figure(figsize=(10,6))
#axs = []
#axs.append(pp.subplot2grid((2,3),(0,0),colspan = 3))
#axs.append(pp.subplot2grid((2,3),(1,0)))
#axs.append(pp.subplot2grid((2,3),(1,1)))
#axs.append(pp.subplot2grid((2,3),(1,2)))

##axs[0].plot(Transient_Times,np.mean(output[:,Node,:].T,axis=1),'k--') 
##axs[0].fill_between(Transient_Times,np.percentile(output[:,Node,:].T,5,axis=1),np.percentile(output[:,Node,:].T,95,axis=1),alpha = 0.5)

##axs[0].plot(Transient_Times,cp.E(foo_approx,demand_distribution),'k')
##axs[0].fill_between(Transient_Times,cp.Perc(foo_approx,5,demand_distribution)[0],cp.Perc(foo_approx,95,demand_distribution)[0],alpha = 0.5)
#axs[0].plot(Transient_Times,np.mean(output[:,Node,:].T,axis=1)-cp.E(foo_approx,distribution),'k') 

##axs[0].plot(Transient_Times,np.percentile(output[:,Node,:].T,5,axis=1)-cp.Perc(foo_approx,5,demand_distribution)[0],'k:') 
##axs[0].plot(Transient_Times,np.percentile(output[:,Node,:].T,95,axis=1)-cp.Perc(foo_approx,95,demand_distribution)[0],'k--') 

#### Initial Time
#a = axs[1].hist(output[:,Node,0],density=True,bins=100,alpha=0.50,color='k')
#x = np.linspace(a[1][0],a[1][-1],1000)
#axs[1].plot(x,coefs_kernal[0].pdf(x),color='k')
#axs[1].set_title('t = 0 s')

#### Middle Time
#a = axs[2].hist(output[:,Node,1500],density=True,bins=100,alpha=0.50,color='k')
#x = np.linspace(a[1][0],a[1][-1],1000)
#axs[2].plot(x,coefs_kernal[1500].pdf(x),color='k')
#axs[2].set_title('t = 15 s')

#### Final Time
#a = axs[3].hist(output[:,Node,-1],density=True,bins=100,alpha=0.50,color='k')
#x = np.linspace(a[1][0],a[1][-1],1000)
#axs[3].plot(x,coefs_kernal[-1].pdf(x),color='k')
#axs[3].set_title('t = 30 s')

#axs[0].set_title('Order: '+str(Order) + ', Collocation Points: ' + str(ConditionPoints))

#axs[0].set_xlim(0,30)
#axs[0].set_xlabel('Time (s)')
#axs[0].set_ylabel('Difference in Expected Head (m)')

#axs[1].set_ylabel('Frequency')
#axs[1].set_xlabel('Head (m)')
#axs[2].set_xlabel('Head (m)')
#axs[3].set_xlabel('Head (m)')
#pp.tight_layout()
#pp.show()

pp.figure(figsize=(10,6))

axs = []
axs.append(pp.subplot2grid((1,2),(0,0)))
axs.append(pp.subplot2grid((1,2),(0,1)))
#axs.append(pp.subplot2grid((2,2),(1,0)))
#axs.append(pp.subplot2grid((2,2),(1,1)))

axs[0].plot(Transient_Times,np.mean(output[:,0,:].T,axis=1),'k') 
axs[0].fill_between(Transient_Times,np.percentile(output[:,0,:].T,5,axis=1),np.percentile(output[:,0,:].T,95,axis=1),color='k',alpha = 0.5)

axs[1].plot(Transient_Times,np.mean(output[:,4,:].T,axis=1),'k') 
axs[1].fill_between(Transient_Times,np.percentile(output[:,4,:].T,5,axis=1),np.percentile(output[:,4,:].T,95,axis=1),color='k',alpha = 0.5)

#axs[2].plot(Transient_Times,cp.E(foo_approx,demand_distribution),'k')
#axs[2].fill_between(Transient_Times,cp.Perc(foo_approx,5,demand_distribution)[0],cp.Perc(foo_approx,95,demand_distribution)[0],color='k',alpha = 0.5)

axs[0].set_xlim(1.5,10)
axs[1].set_xlim(1.5,10)

axs[0].set_xlabel('Time (s)')
axs[1].set_xlabel('Time (s)')
axs[0].set_ylabel('Head (m)')
pp.show()

