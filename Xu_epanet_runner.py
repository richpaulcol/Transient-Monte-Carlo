import numpy as np
import pylab as pp
import csv
from epanettools.epanettools import EPANetSimulation, Node, Link, Network, Nodes,Links, Patterns, Pattern, Controls, Control
from epanettools import epanet2 as et
Directory = 'Projects/Xu_Leak_Detection/'
# FileName = '250701 K709vs2-Export.inp'
FileName = 'ExampleNetwork1.inp'
FileName = 'LeakTestNet.inp'

ret=et.ENopen(Directory+FileName,Directory + 'Temp.rpt',"")
et.ENopenH()
et.ENinitH(0)


time=[]
nodes={}
ret,nnodes=et.ENgetcount(et.EN_NODECOUNT)
for index in range(1,nnodes+1):
    ret,t=et.ENgetnodeid(index)
    nodes[t] = []

for index in range(1,nnodes):
    ret,d = et.ENgetnodevalue(index,et.EN_BASEDEMAND)

    if d == 0.0:
        Pete = np.random.randint(1, 6)
    #if Pete == 1:
        ret1 = et.ENsetnodevalue(index, et.EN_PATTERN, 1)
        ret2 = et.ENsetnodevalue(index,et.EN_BASEDEMAND,300.0)


#Net = EPANetSimulation(Directory+FileName,pdd=False)
#h=Node.value_type['EN_HEAD']
LeakNodes = ['CA81W6J9C0','CBEA1QV88R','CAFPFW92WV','CADXK24N13']
LeakNo = 3

LeakStartTime = 900*int(np.random.uniform(288,960))
print 'Leak Start',LeakStartTime
t=0
while True:
    ### This adds a small amount of variability to the demands
    for index in range(1, nnodes):
        ret, d = et.ENgetnodevalue(index, et.EN_BASEDEMAND)
        ret, dd = et.ENgetnodevalue(index,et.EN_DEMAND)
        newd = float(d)*np.random.uniform(0.99,1.01)
        #newd = np.random.uniform(0,1)
        #if index == 50:
        #print d, newd,dd,index
        ret = et.ENsetnodevalue(index,et.EN_BASEDEMAND,newd)

    ### This adds a small amount of variability to the reservoir
    index = 164
    ret, e = et.ENgetnodevalue(index, et.EN_ELEVATION)
    newe = float(e) * np.random.uniform(0.9995, 1.0005)
    ret = et.ENsetnodevalue(index, et.EN_ELEVATION,newe)


    ret, index = et.ENgetnodeindex(LeakNodes[LeakNo])
    ret = et.ENsetnodevalue(index, et.EN_BASEDEMAND, 0)
    if t == LeakStartTime:
        #print t

        ret = et.ENsetnodevalue(index,et.EN_EMITTER,0.1)
        #print ret


    ret1, t = et.ENrunH()
    time.append(t)
    # for index in range(1, nnodes+1):
    #     ret, node = et.ENgetnodeid(index)
    #     ret, h = et.ENgetnodevalue(index, et.EN_HEAD)
    #     nodes[node].append(h)

    for node in nodes:
        ret, index = et.ENgetnodeindex(node)
        ret, h = et.ENgetnodevalue(index, et.EN_HEAD)
        ret, p = et.ENgetnodevalue(index, et.EN_PRESSURE)
        nodes[node].append(p+np.random.normal(0,0.2))


    ret, tstep = et.ENnextH()


    # ret, index = et.ENgetnodeindex('C9QPT8LCD6')
    # ret, h = et.ENgetnodevalue(index, et.EN_HEAD)
    #print t,h
    if (tstep <= 0):
        break
ret=et.ENcloseH()

sensorList = ['C92ICY6VWN','CBDWGRJVVL','CAKIJ91SJ7','C9EQAVY4YU']

for node in sensorList:
    pp.plot(time,nodes[node])
pp.plot([LeakStartTime,LeakStartTime],[00,100],'k:')
pp.grid('on')
pp.show()


f = open(Directory + 'SensorData_Leak_'+str(LeakNo)+'.csv','wb')
writer = csv.writer(f)
Data = np.array(['Node']+time)
for node in sensorList:
    Data = np.vstack((Data, np.array([node]+nodes[node])))

#np.savetxt(Directory + 'SensorData_Leak_'+str(LeakNo)+'.csv',Data.T,delimiter = ',',fmt='%s %.18e %.18e %.18e %.18e')
writer.writerows(Data.T)
f.close()

f = open(Directory + 'SensorData_' + str(LeakNo) + '_Results.csv','wb')
writer = csv.writer(f)
f.write('The leak position is at '+ LeakNodes[LeakNo] + 'The leak starts at time '+str(LeakStartTime))
f.close()
