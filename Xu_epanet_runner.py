import numpy as np
import pylab as pp
from epanettools.epanettools import EPANetSimulation, Node, Link, Network, Nodes,Links, Patterns, Pattern, Controls, Control
from epanettools import epanet2 as et
Directory = 'Projects/Xu_Leak_Detection/'
FileName = '250701 K709vs2-Export.inp'
FileName = 'ExampleNetwork1.inp'
time=[]
nodes={}
ret,nnodes=et.ENgetcount(et.EN_NODECOUNT)
for index in range(1,nnodes+1):
    ret,t=et.ENgetnodeid(index)
    nodes[t] = []
#Net = EPANetSimulation(Directory+FileName,pdd=False)
#h=Node.value_type['EN_HEAD']
ret=et.ENopen(Directory+FileName,Directory + 'Temp.rpt',"")
et.ENopenH()
et.ENinitH(0)

LeakStartTime = 900*int(np.random.uniform(288,960))
print 'Leak Start',LeakStartTime
while True:
    for index in range(1, nnodes + 1):
        ret, d = et.ENgetnodevalue(index, et.EN_BASEDEMAND)
        newd = float(d)*np.random.uniform(0.8,1.2)
        #print d, newd
        ret = et.ENsetnodevalue(index,et.EN_BASEDEMAND,newd)

    ret, index = et.ENgetnodeindex('C947YS7GPL')
    ret = et.ENsetnodevalue(index, et.EN_BASEDEMAND, 0)
    if t == LeakStartTime:
        print t

        ret = et.ENsetnodevalue(index,et.EN_EMITTER,0.05)
        print ret


    ret1, t = et.ENrunH()
    time.append(t)
    for index in range(1, nnodes+1):
        ret, node = et.ENgetnodeid(index)
        ret, h = et.ENgetnodevalue(index, et.EN_HEAD)
        nodes[node].append(h)
    ret, tstep = et.ENnextH()

    ret, index = et.ENgetnodeindex('C9QPT8LCD6')
    ret, h = et.ENgetnodevalue(index, et.EN_HEAD)
    #print t,h
    if (tstep <= 0):
        break
ret=et.ENcloseH()

for node in nodes:
    pp.plot(time,nodes[node])
pp.show()