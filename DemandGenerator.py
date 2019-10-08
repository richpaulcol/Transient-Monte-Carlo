#import numpy as np
#import pylab as pp
from numpy import zeros, savetxt

Node = 5
Demand = zeros(5000)
Demand[0:500] = 1.
Demand[500:1000] = 2.
Demand[1000:1500] = 3.
Demand[1500:2000] = 4.
Demand[2000:2500] = 5.
Demand[2500:3000] = 4.
Demand[3000:3500] = 3.
Demand[3500:4000] = 2.
Demand[4000:4500] = 1.
Demand[4500:5000] = 0.

savetxt('Demands.csv',Demand,delimiter = ',')
