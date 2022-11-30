import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from model import get_data,spread,stat,statA,get_out,data,T,zero,p,inc_period,asym_period,a,d

data_idx = 0
G,status = get_data(data[data_idx])

history = dict(zip(G.keys(), [[-1,-1,-1] for i in range(1,len(G)+1)]))
origin = np.random.randint(1,len(G),zero)

for i in origin:
	r = np.random.uniform(0,1)
	if r < p:
		rr = np.random.uniform(0,1)
		rr = np.where(asym_period>=rr)[0][0]+1
		status[i] = 'A'+str(rr)
	else:
		rr = np.random.uniform(0,1)
		rr = np.where(inc_period>=rr)[0][0]+1
		status[i] = 'E'+str(rr)
	history[i][0] = 0

for j in range(1,int(T)):
	status,history = spread(G,status,False,a,d,history,j)


fileobj = open('./spreading_data/history','w')
for item in history:
	fileobj.write(str(item) + '\t' +  str(history[item][0])+ '\t' +  str(history[item][1])+ '\t' +  str(history[item][2])+'\n')
fileobj.close()

fileobj = open('./spreading_data/status','w')
for i in status:
	fileobj.write(str(i) + '\t' +  str(status[i])+'\n')
fileobj.close()

get_out(data[0],G,status)
