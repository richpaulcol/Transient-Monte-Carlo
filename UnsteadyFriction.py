import numpy as np
import pylab as pp 
from scipy.interpolate import interp1d


NoNodes = 101
Length = 100
a0 = 500.		#Wave speed
asigma = a0/200.
g = 9.81		#gravity
d0 = 0.05		#mean diameter
dsigma = d0/200.	#std dev diameter
#f = 0.2			#friction
f0 = 0.02
fsigma = f0/200.
Q0 = 0.0005		#initial flow
H0 = 10.		#upstream boundary
maxT = 60.0		#max time in seconds
iters = 1
epsilon = 1.0

#Unsteady Friction Coeffs
#Kut = 0.
#Kux = 1.
k=0.01

xPos = np.linspace(0,Length,NoNodes)
dX = xPos[1]-xPos[0]

Friction  = 'Collins'
Friction  = 'Linear'
#Friction  = 'Collins'
#Friction  = '1Coeff'

hList = []
#pp.figure()
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
	StressList = []
	for t in time:
		
		HNew = np.zeros(len(H))
		QNew = np.zeros(len(Q))
		
		
		### Boundary Conditions
		if t <= 0.5:
			QBC = Q0
		if t>0.5 and t<5.0:
			QBC = QBC + 0.0005 * dT / 4.5
			#QBC = 0.001
		if t>= 5.0:
	#		if QBC > 0:
	#			QBC = 0.0005
	#			#QBC = QBC - 0.00001#0.0005
	#		else:
	#			QBC = 0
				
			QBC = QBC
#		if t > 5:
#			QBC = 0.001
		
		#QBC = Q0
		
		if Friction == 'Linear':
		
			aa = a0
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
			HNew[-1] = (Cp - Bp*Q[-1])
			
			HNew[0] = H0			
			QNew[0] = (H[0] - Cm) / Bm
			
			ShearStress = (1000*d0/4)*(f0*QNew*np.abs(QNew) / (2*g*d0*pipeA[0]**2))
			
#			HNew[ii] = H0
#			QNew[ii] = -(Cm-H0)/Cmm
		
#		elif Friction == 'Chaudry':
#			A = d0**2* np.pi / 4
#			signQdiffQ = np.sign(np.gradient(Q)*Q)
#			###	Unsteady Friction Relationships
#			alpha_p = 2*(1+Kut) / (-signQdiffQ[:-1]*Kux + 2 + Kut)
#			alpha_m = 2*(1+Kut) / (-signQdiffQ[1:]*Kux - 2 - Kut)
#	#		alpha_p = 2*a / (-signQdiffQ[:-1]*Kux + 2 + Kut)
#	#		alpha_m = 2*a / (-signQdiffQ[1:]*Kux - 2 - Kut)
#	
#			for ii in range(H.size):
#				if ii == 0: ###Upstream Boundary Condition
##					if Q[ii]*(Q[ii+1]-Q[ii]) >=0:
##						S = 1
##					else:
##						S = -1
##						
#					
#					
#					S = np.sign(Q[ii]*(Q[ii+1]-Q[ii]))
#					alpha_m = 2*a0*(1+Kut)/(-S*Kux- 2 - Kut)
#					
#					Hm = H[ii] - alpha_m * dT/dX *( H[ii+1] - H[ii])
#					Qm = Q[ii] - alpha_m * dT/dX *( Q[ii+1] - Q[ii])
#					
#					Ca = g*A/a0/alpha_m
#					Cn = Qm + Ca * Hm - f0 * Qm*abs(Qm) * dT / (2*d0*(1+Kut))
#					
#					HNew[ii] = H0
#					QNew[ii] = Cn + Ca * H0
#				
#				elif ii == H.size-1: ### Downstream Boundary Condition
#					
##					if Q[ii]*(Q[ii]-Q[ii-1]) >=0:
##						S = 1
##					else:
##						S = -1
#					S = np.sign(Q[ii]*(Q[ii]-Q[ii-1]))
#					alpha_p = 2*a0*(1+Kut)/(-S*Kux + 2 + Kut)
#					
#					Hp = H[ii] - alpha_p * dT/dX *( H[ii] - H[ii-1])
#					Qp = Q[ii] - alpha_p* dT/dX *( Q[ii] - Q[ii-1])
#					
#					Ca = g*A/a0/alpha_p
#					Cp = Qp + Ca * Hp - f0 * Qp*abs(Qp) * dT / (2*d0*(1+Kut))
#					QNew[ii] = QBC
#					HNew[ii] = (QBC - Cp)/(-Ca)
#					#print alpha_p
#				else:
##					if Q[ii]*(Q[ii+1]-Q[ii-1]) >=0:
##						S = 1
##					else:
##						S =-1
#					S = np.sign(Q[ii]*(Q[ii+1]-Q[ii-1]))
#					alpha_m = 2*a0*(1+Kut)/(-S*Kux- 2 - Kut)
#					alpha_p = 2*a0*(1+Kut)/(-S*Kux + 2 + Kut)
#					
#					Hm = H[ii] - alpha_m * dT/dX *( H[ii+1] - H[ii])
#					Qm = Q[ii] - alpha_m * dT/dX *( Q[ii+1] - Q[ii])
#					Hp = H[ii] - alpha_p * dT/dX *( H[ii] - H[ii-1])
#					Qp = Q[ii] - alpha_p* dT/dX *( Q[ii] - Q[ii-1])
#					
#					Can = g*A/a0/alpha_m
#					Cn = Qm + Can * Hm - f0 * Qm*abs(Qm) * dT / (2*d0*(1+Kut))
#					Cap = g*A/a0/alpha_p
#					Cp = Qp + Cap * Hp - f0 * Qp*abs(Qp) * dT / (2*d0*(1+Kut))
#					
#					QNew[ii] = (Cp*Can + Cn*Cap) /  (Cap+Can)
#					HNew[ii] = (-1./2)*(Cp - Cn) / (Cap+Can)
#			if t == 0:
#				asfasfsa
				
				
		elif Friction == '1Coeff':
			A = d0**2* np.pi / 4
			for ii in range(H.size):
				if ii == 0: ###Upstream Boundary Condition
					if Q[ii]*(Q[ii+1]-Q[ii]) >=0:
						S = 1
					else:
						S = -1
						
					
					S = np.sign(Q[ii]*(Q[ii+1]-Q[ii]))
					beta_m = (1.+k)/(1.+(1./2.)*k*(1.+S))
					
					Hb = H[ii] - beta_m *a0/(1+k)*dT/dX *( H[ii] - H[ii+1])
					Qb = Q[ii] - beta_m *a0/(1+k)*dT/dX *( Q[ii] - Q[ii+1])
					
					Cmm = ((1+k)/beta_m)*abs(Qb)*dX*f0/(2*g*d0*A**2)+beta_m*a0/(A*g)
					Cm = Hb-beta_m*a0*Qb/(A*g)
					
					HNew[ii] = H0
					QNew[ii] = -(Cm-H0)/Cmm
				
				elif ii == H.size-1: ### Downstream Boundary Condition
					
					if Q[ii]*(Q[ii]-Q[ii-1]) >=0:
						S = 1
					else:
						S = -1
					S = np.sign(Q[ii]*(Q[ii]-Q[ii-1]))
					beta_p = (1.+k)/(1.+(1/2.)*k*(1.-S))
					Ha = H[ii] - beta_p *a0/(1+k)* dT/dX *( H[ii] - H[ii-1])
					Qa = Q[ii] - beta_p *a0/(1+k)* dT/dX *( Q[ii] - Q[ii-1])
					
					Cpp = -((1+k)/beta_p)*abs(Qa)*dX*f0/(2*g*d0*A**2)-beta_p*a0/(A*g)
					Cp = Ha+beta_p*a0*Qa/(A*g)
					
					QNew[ii] = QBC
					HNew[ii] = Cpp*QBC+Cp
					#print alpha_p
				else:
					if Q[ii]*(Q[ii+1]-Q[ii-1]) >=0:
						S = 1
					else:
						S =-1
					S = np.sign(Q[ii]*(Q[ii+1]-Q[ii-1]))
					beta_m = (1.+k)/(1.+(1./2.)*k*(1.+S))
					beta_p = (1.+k)/(1.+(1/2.)*k*(1.-S))
					
					Hb = H[ii] - beta_m *a0/(1+k)*dT/dX *( H[ii] - H[ii+1])
					Qb = Q[ii] - beta_m *a0/(1+k)* dT/dX *( Q[ii] - Q[ii+1])
					Ha = H[ii] - beta_p *a0/(1+k)* dT/dX *( H[ii] - H[ii-1])
					Qa = Q[ii] - beta_p *a0/(1+k)* dT/dX *( Q[ii] - Q[ii-1])
					
					Cmm = ((1+k)/beta_m)*abs(Qb)*dX*f0/(2*g*d0*A**2)+beta_m*a0/(A*g)
					Cm = Hb-beta_m*a0*Qb/(A*g)
					Cpp = -((1+k)/beta_p)*abs(Qa)*dX*f0/(2*g*d0*A**2)-beta_p*a0/(A*g)
					Cp = Ha+beta_p*a0*Qa/(A*g)
					
					QNew[ii] = -(Cm-Cp)/(Cmm-Cpp)
					HNew[ii] = -(Cm*Cpp-Cmm*Cp)/(Cmm-Cpp)
#			if t ==0:
#				Dave	
			ShearStress = (1000*d0/4)*((f0*QNew*np.abs(QNew) / (2*g*d0*pipeA[0]**2)) + (k/g/A)*((QNew-Q)/dT + a0 * np.gradient(Q)/dX))
				
#		elif Friction == 'Collins':
#			A = d0**2* np.pi / 4
#			S = np.sign(np.gradient(Q)*Q)
#			R = f*dX / (2 * g * pipeD*pipeA**2)
#			beta_p = (-Kux * S + np.sqrt(Kux**2 + 4 + 4.*Kut)) / (Kut +1)
#			beta_m = (-Kux * S - np.sqrt(Kux**2 + 4 + 4.*Kut)) / (Kut +1)
#			
#			Cpp = g*A/(2*a0) * beta_p
#			Cmm = g*A/(2*a0) * beta_m
#			
#			aa = a0
#			Hp = H[1:-1] - aa * dT/dX *( H[1:-1] - H[:-2])
#			Hm = H[1:-1] - aa * dT/dX *( H[1:-1] - H[2:])
#			
#			Qp = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[:-2])
#			Qm = Q[1:-1] - aa * dT/dX *( Q[1:-1] - Q[2:])
#			
#			Cp = Qp + Cpp[0:-2] * Hp - R[0:-1]*dT*Qp*abs(Qp)
#			Cm = Qm + Cmm[2:] * Hm - R[1:]*dT*Qm*abs(Qm)
#			
#			for ii in range(H.size):
#				if ii == 0: ###Upstream Boundary Condition
#				
#					
#				
#					if Q[ii]*(Q[ii+1]-Q[ii]) >=0:
#						S = 1
#					else:
#						S = -1
#						
#					
#					
#					#S = np.sign(Q[ii]*(Q[ii+1]-Q[ii]))
#					alpha_m = 2*a0*(1+Kut)/(-S*Kux- 2 - Kut)
#					
#					Hm = H[ii] - alpha_m * dT/dX *( H[ii+1] - H[ii])
#					Qm = Q[ii] - alpha_m * dT/dX *( Q[ii+1] - Q[ii])
#					
#					Ca = g*A/a0/alpha_m
#					Cn = Qm + Ca * Hm - f0 * Qm*abs(Qm) * dT / (2*d0*(1+Kut))
#					
#					HNew[ii] = H0
#					QNew[ii] = Cn + Ca * H0
#				
#				elif ii == H.size-1: ### Downstream Boundary Condition
#					
#					if Q[ii]*(Q[ii]-Q[ii-1]) >=0:
#						S = 1
#					else:
#						S = -1
##					S = np.sign(Q[ii]*(Q[ii]-Q[ii-1]))
#					alpha_p = 2*a0*(1+Kut)/(-S*Kux + 2 + Kut)
#					
#					Hp = H[ii] - alpha_p * dT/dX *( H[ii] - H[ii-1])
#					Qp = Q[ii] - alpha_p* dT/dX *( Q[ii] - Q[ii-1])
#					
#					Ca = g*A/a0/alpha_p
#					Cp = Qp + Ca * Hp - f0 * Qp*abs(Qp) * dT / (2*d0*(1+Kut))
#					QNew[ii] = QBC
#					HNew[ii] = (QBC - Cp)/(-Ca)
#					#print alpha_p
#				else:
#					if Q[ii]*(Q[ii+1]-Q[ii-1]) >=0:
#						S = 1
#					else:
#						S =-1
##					S = np.sign(Q[ii]*(Q[ii+1]-Q[ii-1]))
#					alpha_m = 2*a0*(1+Kut)/(-S*Kux- 2 - Kut)
#					alpha_p = 2*a0*(1+Kut)/(-S*Kux + 2 + Kut)
#					
#					Hm = H[ii] - alpha_m * dT/dX *( H[ii+1] - H[ii])
#					Qm = Q[ii] - alpha_m * dT/dX *( Q[ii+1] - Q[ii])
#					Hp = H[ii] - alpha_p * dT/dX *( H[ii] - H[ii-1])
#					Qp = Q[ii] - alpha_p* dT/dX *( Q[ii] - Q[ii-1])
#					
#					Can = g*A/a0/alpha_m
#					Cn = Qm + Can * Hm - f0 * Qm*abs(Qm) * dT / (2*d0*(1+Kut))
#					Cap = g*A/a0/alpha_p
#					Cp = Qp + Cap * Hp - f0 * Qp*abs(Qp) * dT / (2*d0*(1+Kut))
#					
#					QNew[ii] = (Cp*Can + Cn*Cap) /  (Cap+Can)
#					HNew[ii] = (-1./2)*(Cp - Cn) / (Cap+Can)
#			if t == 0:
#				asfasfsa
		
		H = np.copy(HNew)
		Q = np.copy(QNew)
		
		
		QList.append(Q)
		HList.append(H)
		StressList.append(ShearStress)
		#break
		
	QArray = np.array(QList)
	HArray = np.array(HList)
	QData.append(QArray)
	HData.append(HArray)
	StressArray = np.array(StressList)
	#hList.append(H)
#	pp.plot(time,HArray[:,0])
#	pp.plot(time,HArray[:,NoNodes/10])
	f,axs = pp.subplots(3,1,sharex=True)
	axs[0].plot(time,HArray[:,NoNodes/2])
	axs[1].plot(time,QArray[:,NoNodes/2])
	axs[2].plot(time,StressArray[:,NoNodes/2])
	
	if Friction == 'Linear':
		Directory = 'Projects/UncertainShear/'
#		np.save(Directory+'LinearHeadFast.npy',HArray)
#		np.save(Directory+'LinearFlowFast.npy',QArray)
#		np.save(Directory+'LinearShearFast.npy',StressArray)
		
		np.save(Directory+'LinearHeadSlow.npy',HArray)
		np.save(Directory+'LinearFlowSlow.npy',QArray)
		np.save(Directory+'LinearShearSlow.npy',StressArray)
		
	if Friction == '1Coeff':
		Directory = 'Projects/UncertainShear/'
#		np.save(Directory+'USHeadFast.npy',HArray)
#		np.save(Directory+'USFlowFast.npy',QArray)
#		np.save(Directory+'USShearFast.npy',StressArray)
		
		np.save(Directory+'USHeadSlow.npy',HArray)
		np.save(Directory+'USFlowSlow.npy',QArray)
		np.save(Directory+'USShearSlow.npy',StressArray)
#	pp.plot(time,HArray[:,-1])
pp.show()	
#QData = np.array(QData)
HData = np.array(HData)

#np.save('Flow SectionL ' + str(SectionL) + ' Iters ' + str(iters),QData)
#np.save('Head SectionL ' + str(SectionL) + ' Iters ' + str(iters),HData)




pp.show()
