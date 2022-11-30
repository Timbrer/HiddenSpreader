import numpy as np
import copy
import scipy.stats as stats    
import networkx as nx
import matplotlib.pyplot as plt
from model import get_data,get_status,stat,statA,data,a,d,T as ori_T,p

data_idx = 0
beta = a+d
T = ori_T

inc_prob = 1-stats.lognorm.cdf(np.arange(0,999)/4.17,s=0.66)
remove_prob = 1-stats.norm.cdf((np.arange(0,999)-8.83)/3.88)
asym_prob = 1-stats.norm.cdf((np.arange(0,999)-20)/5)

inc_prob_ra = [inc_prob[i]/inc_prob[i-1] for i in range(1,40)]
remove_prob_ra = [remove_prob[i]/remove_prob[i-1] for i in range(1,40)]
asym_prob_ra = [asym_prob[i]/asym_prob[i-1] for i in range(1,40)]

def load_data(status):
	idx_features_labels1 = np.genfromtxt(r'./spreading_data/train', dtype=np.int32)
	try:
		idx_train = set(idx_features_labels1[:,0])
	except:
		idx_train = set([idx_features_labels1[0]])
	
	return idx_train
	
def simulate_spread(state,G,idx_train,nodes):
	
	new_state = {}
	for i in nodes:
		new_state[i] = np.empty(77)
	for i in idx_train:
		new_state[i] = state[i]
	''' predict: 0S 1-20E 21-40I 41-75A 76R'''
	for t in range(1,T):
		for i in nodes:
			new_state[i][-1] = state[i][-1] + state[i][40] + state[i][75]
			new_state[i][21] = state[i][20]
			for j in range(2,21):
				new_state[i][j] = state[i][j-1]*inc_prob_ra[j-2]
				new_state[i][21] += state[i][j-1]*(1-inc_prob_ra[j-2])
			for j in range(22,41):
				new_state[i][j] = state[i][j-1]*remove_prob_ra[j-22]
				new_state[i][-1] += state[i][j-1]*(1-remove_prob_ra[j-22])
			for j in range(42,76):
				new_state[i][j] = state[i][j-1]*asym_prob_ra[j-42]
				new_state[i][-1] += state[i][j-1]*(1-asym_prob_ra[j-42])
			
			neighbor_inf = 1
			for nei in G[i]:
				if len(state[nei]) == 2:
					if t>state[nei][0] and t <= state[nei][1]:
						neighbor_inf *= (1-beta)
				else:
					neighbor_inf *= (1-beta*(sum(state[nei][1:21])+sum(state[nei][41:76])))
			
			new_state[i][1] = state[i][0]*(1-p)*(1-neighbor_inf)
			new_state[i][41] = state[i][0]*p*(1-neighbor_inf)
			new_state[i][0] = state[i][0]*neighbor_inf

		state = copy.deepcopy(new_state)
	
	return state
	
	
G,status = get_data(data[data_idx])
idx_train = load_data(status)

nodes = set(list(G.keys()))-idx_train

status = np.genfromtxt(r"./spreading_data/status", dtype=np.dtype(str))
status = dict(zip(status[:,0].astype(np.int32), status[:,-1]))

history = np.genfromtxt(r"./spreading_data/history", dtype=np.int32)
history = dict(zip(history[:,0],history[:,1:]))

state = {}
for i in G:
	if i in idx_train:
		if history[i][1] == -1:
			state[i] = (history[i][0],T)
		else:
			state[i] = (history[i][0],history[i][1])
	else:
		state[i] = np.zeros(77)
		state[i][0] = 1

state = simulate_spread(state,G,idx_train,nodes)


ans = np.array([ [i,sum(state[i][1:21])+sum(state[i][41:76])] for i in nodes ])

temp_ans = np.argsort(ans[:,1])

fileobj = open(r'./out.txt','w')

for i in temp_ans:
	fileobj.write(str(int(ans[i][0]))+'\t'+str(ans[i][1])+'\t'+str(int(status[ans[i][0]]!='S'))+'\n')
			
fileobj.close()

edges = np.genfromtxt("./data/"+data[data_idx], dtype=np.int32)[:,0:2]
G=nx.Graph()
G.add_edges_from(list(edges))
pos = nx.spring_layout(G)
sym_patient = np.array([[i,np.max(ans[:,1])] for i in idx_train])
all_node_color = np.vstack((ans,sym_patient)) 
nx.draw(G, pos, node_color=all_node_color[np.argsort(all_node_color[:,0])][:,1], node_size=200, cmap=plt.cm.RdYlGn)
plt.savefig('rank.svg')
