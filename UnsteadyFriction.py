import numpy as np
import pylab as pp 

import numpy as np
import pylab as pp

NoNodes = 101
Length = 100
a0 = 250.		#Wave speed
asigma = a0/200.
g = 9.81		#gravity
d0 = 0.05		#mean diameter
dsigma = d0/200.	#std dev diameter
f = 0.2			#friction
f0 = 0.2
fsigma = f0/200.
Q0 = 0.001		#initial flow
H0 = 10.		#upstream boundary
maxT = 10.0		#max time in seconds
iters = 1
epsilon = 0.85

xPos = np.linspace(0,Length,NoNodes)
dX = xPos[1]-xPos[0]


hList = []
pp.figure()
QData = []
HData = []
for i in range(iters):
	print i
	#### setting out the size of the arrays for the pipe properties
	pipeD = np.zeros(len(xPos)-1)	## Looking at pipe diameters
	f = np.zeros(len(xPos)-1)
	a = np.zeros(len(xPos)-1)
	
	### adding a random element to the pipe properties
	#for i in range(0,len(pipeD)/SectionL):
		#pipeD[i*SectionL:i*SectionL+SectionL] = pp.normal(d0,dsigma)
		#f[i*SectionL:i*SectionL+SectionL] = pp.normal(f0,fsigma)
		#a[i*SectionL:i*SectionL+SectionL] =pp.normal(a0,asigma)
	
	###  Populating the arrays with the given values for this simulation
	a = np.ones(a.shape)*a0
	pipeD = np.ones(pipeD.shape)*d0
	pipeA = pipeD**2 * np.pi / 4.	
	
	###  Setting the time step
	dT = dX / a0#pp.mean(a)
	#print pipeD
	
	### Initating Head and Flow arrays
	H = np.zeros(len(xPos))
	Q = np.ones(len(xPos))*Q0	## Initial Flow conditon
	
	#Setting the initial conditions of the pressure and 
	H[0] = H0
	for i in range(1,len(H)):
		H[i] = H[i-1] - f[i-1] *dX* Q[i-1]**2 / (2. * g *pipeA[i-1]**2 *pipeD[i-1])
	Hinitial = H
	HRight = H[-1]
	time = np.arange(0,maxT+dT,dT)	
	
	B = a / (g * pipeA)
	R = f*dX / (2 * g * pipeD*pipeA**2)
	
	QList = []
	HList = []
	for t in time:
		
		HNew = np.zeros(len(H))
		QNew = np.zeros(len(Q))
		
		Cp = H[:-2] + Q[:-2]*(B[:-1] - R[:-1]*abs(Q[:-2])*(1-epsilon))
		Cm = H[2:] - Q[2:]*(B[1:] - R[1:]*abs(Q[2:])*(1-epsilon)) 
		
		Bp = B[:-1] + epsilon*R[:-1]*abs(Q[:-2])
		Bm = B[1:] + epsilon*R[1:]*abs(Q[2:])
		
		#QNew[1:-1] = (Cp[1:] - Cm[:-1])/(Bp[1:] + Bm[:-1])
		QNew[1:-1] = (Cp - Cm)/(Bp + Bm)
		#HNew[1:-1] = (Cm[:-1]*Bp[1:] + Bm[:-1]*Cp[1:])/(Bp[1:]+Bm[:-1])
		HNew[1:-1] = (Cm*Bp + Bm*Cp)/(Bp+Bm)
		
		#if t > 0.5:
		#	Q0 = 0.
		
		if t > 0.5:
			if QBC > 0:
				QBC = 0.0005
				#QBC = QBC - 0.00001#0.0005
			else:
				QBC = 0
		if t > 5:
			QBC = 0.001
		if t <= 0.5:
			QBC = Q0
		
		Cp = H[-2] + Q[-2]*(B[-1] - R[-1]*abs(Q[-2])*(1-epsilon))
		Cm = H[1] - Q[1]*(B[0] - R[0]*abs(Q[1])*(1-epsilon))
		Bp = B[-1] + epsilon*R[-1]*abs(Q[-2])
		Bm = B[0] + epsilon*R[0]*abs(Q[1])		
		
		QNew[-1] = QBC#(Cp - H[-1]) / Bp
		HNew[0] = H0
		HNew[-1] = (Cp - Bp*Q[-1])
		QNew[0] = (H[0] - Cm) / Bm
		
		H = np.copy(HNew)
		Q = np.copy(QNew)
		
		QList.append(Q)
		HList.append(H)
		#break
		
	QArray = np.array(QList)
	HArray = np.array(HList)
	QData.append(QArray)
	HData.append(HArray)
	#hList.append(H)
	pp.plot(time,HArray[:,10])
	pp.plot(time,HArray[:,NoNodes/2])
	pp.plot(time,HArray[:,-1])
pp.show()	
#QData = np.array(QData)
HData = np.array(HData)

#np.save('Flow SectionL ' + str(SectionL) + ' Iters ' + str(iters),QData)
#np.save('Head SectionL ' + str(SectionL) + ' Iters ' + str(iters),HData)

pp.show()
