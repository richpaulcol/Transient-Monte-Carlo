import numpy as np
import pylab as pp

dH = 100
Lambda = 0.02
l = 10000
d = 0.3
A = d**2 * np.pi / 4

dO = np.linspace(0.001,d,100)
AO = dO**2 * np.pi / 4
Cd = 0.66

R = (Lambda * l * (AO**2 / A**2) / d + (1./Cd**2))

VO = np.sqrt(2*9.81*dH / R)

pp.plot(AO/A,VO/np.sqrt(2*9.81*dH),label = 'Upstream Length = ' + str(l) + 'm')
pp.xlabel('Burst Area / Pipe Area')
pp.ylabel('Burst Velocity (m/s) / (2 g dH)^0.5')


#Rb = 1.8 *(Q * H/(g**0.5 * Cd * d0**(7/2)))**0.243

