import numpy as np
import pylab as pp 
from scipy.interpolate import interp1d


NoNodes = 101
Length = 100
a0 = 250.		#Wave speed
asigma = a0/200.
g = 9.81		#gravity
d0 = 0.05		#mean diameter
dsigma = d0/200.	#std dev diameter
#f = 0.2			#friction
f0 = 0.2
fsigma = f0/200.
Q0 = 0.001		#initial flow
H0 = 10.		#upstream boundary
maxT = 10.0		#max time in seconds
iters = 1
epsilon = 0.85

#Unsteady Friction Coeffs
Kut = 0.00
Kux = 0.00

xPos = np.linspace(0,Length,NoNodes)
dX = xPos[1]-xPos[0]

Friction  = 'Collins'


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
	f = np.ones(f.shape)*f0
	pipeA = pipeD**2 * np.pi / 4.
		
	
	###  Setting the time step
	dT = dX / a0#pp.mean(a)
	#print pipeD
	
	### Initating Head and Flow arrays
	H = np.zeros(len(xPos))
	Q = np.ones(len(xPos))*Q0	## Initial Flow conditon
	
	signQdiffQ = np.sign(np.gradient(Q)*Q)	## Initializing the QdiffQ vector
	
	#Setting the initial conditions of the pressure head
	H[0] = H0
	for i in range(1,len(H)):
		H[i] = H[i-1] - f[i-1] *dX* Q[i-1]**2 / (2. * g *pipeA[i-1]**2 *pipeD[i-1])
	#Hinitial = H
	#HRight = H[-1]
	
	### Intitiating the time vector
	time = np.arange(0,maxT+dT,dT)	
	
	
	
	
	QList = []
	HList = []
	for t in time:
		
		HNew = np.zeros(len(H))
		QNew = np.zeros(len(Q))
		
		
		### Boundary Conditions
		if t > 0.5:
	#		if QBC > 0:
	#			QBC = 0.0005
	#			#QBC = QBC - 0.00001#0.0005
	#		else:
	#			QBC = 0
				
			QBC = 0.002
		if t > 5:
			QBC = 0.001
		if t <= 0.5:
			QBC = Q0
		
		
		if Friction == 'Linear':
		
			aa = 250
			Hp = H[1:-1] - aa * dT/dX *( H[1:-1] - H[:-2])
			Hm = H[1:-1] - aa * dT/dX *( H[1:-1] - H[2:])
			
			Qp = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[:-2])
			Qm = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[2:])
			
			B = aa / (g * pipeA)
			R = f*dX / (2 * g * pipeD*pipeA**2)
			
			
			### These are the version with interpolation
			Cp = Hp + Qp*(B[:-1] - R[:-1]*abs(Qp)*(1-epsilon))
			Cm = Hm - Qm*(B[1:] - R[1:]*abs(Qm)*(1-epsilon)) 
			Bp = B[:-1] + epsilon*R[:-1]*abs(Qp)
			Bm = B[1:] + epsilon*R[1:]*abs(Qm)
		
		
		
			### These are the version without interpolation
	#		Cp = H[:-2] + Q[:-2]*(B[:-1] - R[:-1]*abs(Q[:-2])*(1-epsilon))
	#		Cm = H[2:] - Q[2:]*(B[1:] - R[1:]*abs(Q[2:])*(1-epsilon)) 
	#		
	#		Bp = B[:-1] + epsilon*R[:-1]*abs(Q[:-2])
	#		Bm = B[1:] + epsilon*R[1:]*abs(Q[2:])

			QNew[1:-1] = (Cp - Cm)/(Bp + Bm)
			HNew[1:-1] = (Cm*Bp + Bm*Cp)/(Bp+Bm)
			
			Cp = H[-2] + Q[-2]*(B[-1] - R[-1]*abs(Q[-2])*(1-epsilon))
			Cm = H[1] - Q[1]*(B[0] - R[0]*abs(Q[1])*(1-epsilon))
			Bp = B[-1] + epsilon*R[-1]*abs(Q[-2])
			Bm = B[0] + epsilon*R[0]*abs(Q[1])		
		
		
			QNew[-1] = QBC#(Cp - H[-1]) / Bp
			HNew[0] = H0
			HNew[-1] = (Cp - Bp*Q[-1])
			QNew[0] = (H[0] - Cm) / Bm
		
		elif Friction == 'Chaudry':
			signQdiffQ = np.sign(np.gradient(Q)*Q)
			###	Unsteady Friction Relationships
			alpha_p = 2*(1+Kut) / (-signQdiffQ[:-1]*Kux + 2 + Kut)
			alpha_m = 2*(1+Kut) / (-signQdiffQ[1:]*Kux - 2 - Kut)
	#		alpha_p = 2*a / (-signQdiffQ[:-1]*Kux + 2 + Kut)
	#		alpha_m = 2*a / (-signQdiffQ[1:]*Kux - 2 - Kut)
			
			
			Hp = H[1:-1] - a[1:]*alpha_p[1:] * dT/dX *( H[1:-1] - H[:-2])
			Hm = H[1:-1] + a[:-1]*alpha_m[:-1] * dT/dX *( H[1:-1] - H[2:])
			
			Qp = Q[1:-1] - a[1:]*alpha_p[1:] * dT/dX *( Q[1:-1] - Q[:-2])
			Qm = Q[1:-1] + a[:-1]*alpha_m[:-1] * dT/dX *( Q[1:-1] - Q[2:])
			
			Cap = (g*pipeA / a) / alpha_p
			Cpp = Qp + Cap[:-1] * Hp  - f[:-1] * Qp * np.abs(Qp) * dT/ (2*pipeD[:-1]*(1+Kut))
			
			Cam = (g*pipeA / a) / alpha_m
			Cmm = Qm + Cam[1:]*Hm - f[1:] * Qm * np.abs(Qm) *dT/ (2*pipeD[1:]*(1+Kut))
		
			HNew[1:-1] = (Cpp - Cmm)/(Cap[:-1] - Cam[1:])
			QNew[1:-1] = (Cmm*Cap[:-1] - Cpp*Cam[1:])/(Cap[:-1]-Cam[1:])
			
	#		sadasd
			
			
				
				
			##
			
			
			### Unsteady Friction Versions of the Boundary Conditions
			alpha_p = 2*(1+Kut) / (-signQdiffQ[-1]*Kux + 2 + Kut)
			alpha_m = 2*(1+Kut) / (-signQdiffQ[0]*Kux - 2 - Kut)
	#		alpha_p = 2*a[-1] / (-signQdiffQ[-1]*Kux + 2 + Kut)
	#		alpha_m = 2*a[0] / (-signQdiffQ[0]*Kux - 2 - Kut)
			
			
			Hp = H[-1] -a[-1]*alpha_p * dT/dX *( H[-1] - H[-2])
			Hm = H[0] + a[-1]*alpha_m * dT/dX *( H[0] - H[1])
			Qp = Q[-1] - a[0]*alpha_p * dT/dX *( Q[-1] - Q[-2])
			Qm = Q[0] + a[0]*alpha_m * dT/dX *( Q[0] - Q[1])
			
			Cap = (g*pipeA[-1] / a[-1]) / alpha_p
			Cpp = Qp + Cap * Hp  - f[-1] * Qp * np.abs(Qp) *dT/ (2*pipeD[-1]*(1+Kut))
			
			Cam = (g*pipeA[0] / a[0]) / alpha_m
			Cmm = Qm + Cam*Hm - f[0] * Qm * np.abs(Qm) *dT/ (2*pipeD[0]*(1+Kut))
			
			
			HNew[0] = H0
			QNew[0] = Cmm + Cam*H0
					
			QNew[-1] = QBC#(Cp - H[-1]) / Bp
			HNew[-1] = (QNew[-1] - Cpp)/(-Cap)
				
		elif Friction == 'Collins':
			A = d0**2* np.pi / 4
			S = np.sign(np.gradient(Q)*Q)
			R = f*dX / (2 * g * pipeD*pipeA**2)
			beta_p = (-Kux * S + np.sqrt(Kux**2 + 4 + 4.*Kut)) / (Kut +1)
			beta_m = (-Kux * S - np.sqrt(Kux**2 + 4 + 4.*Kut)) / (Kut +1)
			
			Cpp = g*A/(2*a0) * beta_p
			Cmm = g*A/(2*a0) * beta_m
			
			aa = a0
			Hp = H[1:-1] - aa * dT/dX *( H[1:-1] - H[:-2])
			Hm = H[1:-1] - aa * dT/dX *( H[1:-1] - H[2:])
			
			Qp = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[:-2])
			Qm = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[2:])
			
			Cp = Qp + Cpp[0:-2] * Hp - R[0:-1]*dT*Qp*abs(Qp)
			Cm = Qm + Cmm[2:] * Hm - R[1:]*dT*Qm*abs(Qm)
			
			QNew[1:-1] = -(Cm*Cpp[:-2]-Cmm[2:]*Cp)/(Cmm[:-2]-Cpp[2:])
			HNew[1:-1] = -(Cm-Cp)/(Cmm[:-2]-Cpp[2:])
			
			HNew[0] = H0
			QNew[0] = Q[1] + Cmm[0] * H[1] - R[0]*dT*Q[1]*abs(Q[1]) + Cmm[0]*H0
					
			QNew[-1] = QBC#(Cp - H[-1]) / Bp
			HNew[-1] = (QBC - Q[-2] + Cpp[-1] * H[-2] - R[-1]*dT*Q[-2]*abs(Q[-2]))/(-Cpp[-1])
		
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
	pp.plot(time,HArray[:,0])
	pp.plot(time,HArray[:,NoNodes/10])
	pp.plot(time,HArray[:,NoNodes/2])
	pp.plot(time,HArray[:,-1])
pp.show()	
#QData = np.array(QData)
HData = np.array(HData)

#np.save('Flow SectionL ' + str(SectionL) + ' Iters ' + str(iters),QData)
#np.save('Head SectionL ' + str(SectionL) + ' Iters ' + str(iters),HData)

pp.show()
