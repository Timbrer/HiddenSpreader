import numpy as np
import scipy.sparse as sp
import copy
import math
import scipy.stats as stats    

data = ['infectious','univeristy_email']
zero = 5                     # Number of zero infections 
k = 13.49                    # Average degree                                                  
R0 = 3.5                     # Reproduction number
TE = 5.18                    # Mean presymptomatic period 
TI = 8.83                    # Mean symptomatic period 
TA = 20                      # Mean asymptomatic period 
p = 0.15                     # Ratio of asymptomatic patience
T=11                         # Spreading days
E_L_Time =  (1-p) * (TE+TI) + p*TA  # Mean virus-carry period
a = R0/(k*E_L_Time) *(1-p)   
d = R0/(k*E_L_Time) *p

inc_period = stats.lognorm.cdf(np.arange(1,41)/4.17,s=0.66)
inc_period[-1] = 1
remove_period = stats.norm.cdf((np.arange(1,31)-8.82)/3.88)
remove_period[-1] = 1
asym_period = stats.norm.cdf((np.arange(1,41)-20)/5)
asym_period[-1] = 1

def inc(time):
	return stats.lognorm.pdf(x/4.17,s=0.66)/4.17

def remove(time):
	return stats.norm.pdf((x-8.83)/3.88)/3.88

def asym(x):
	return stats.norm.pdf((x-20)/5)/5

def get_status(name):
	status = np.genfromtxt(r"./spreading_data/status", dtype=np.dtype(str))
	return dict(zip(status[:,0].astype(np.int32), status[:,-1]))
	
	
def get_data(name):
	edge = np.genfromtxt(r"./data/{}".format(name), dtype=np.int32)
	G = {}
	for (i,j) in edge[:,0:2]:
		if i in G:
			G[i].add(j)
		else:
			G[i] = set([j])
		if j in G:
			G[j].add(i)
		else:
			G[j] = set([i])
	return G, dict(zip(G.keys(),['S']*len(G)))

def spread(G,status,control,a,d,history,T):
	new_status = copy.deepcopy(status)
	for i in status:
		if status[i][0] == 'E' or status[i][0] == 'I' or status[i][0] == 'A':
			if status[i] != 'IR' and status[i] != 'AR':
				for j in G[i]:
					if status[j] == 'S':
						#print(j)
						r = np.random.uniform(0,1)
						if r < a:
							rr = np.random.uniform(0,1)
							rr = np.where(inc_period>=rr)[0][0]+1
							new_status[j] = 'E'+str(rr)
							if history[j][2] == -1:
								history[j][2] = T-history[i][2]
							elif history[j][2] > T-history[i][2]:
								history[j][2] = T-history[i][2]
							history[j][0] = T
						elif r>a and r<a+d:
							rr = np.random.uniform(0,1)
							rr = np.where(asym_period>=rr)[0][0]+1
							new_status[j] = 'A'+str(rr)
							if history[j][2] == -1:
								history[j][2] = T-history[i][2]
							elif history[j][2] > T-history[i][2]:
								history[j][2] = T-history[i][2]
							history[j][0] = T
				if status[i][0] == 'E':
						if T-history[i][0] == int(status[i][1:]):
							rr = np.random.uniform(0,1)
							rr = np.where(remove_period>=rr)[0][0]+1
							new_status[i] = 'I'+str(rr+int(status[i][1:]))
				elif status[i][0] == 'I':
					if control == False:
						if T-history[i][0] == int(status[i][1:]):
							new_status[i] = 'IR'
							history[i][1] = T
					elif control == True:
						new_status[i] = 'IR'
						history[i][1] = T
				elif status[i][0] == 'A':
					if T-history[i][0] == int(status[i][1:]):
						new_status[i] = 'AR'
						history[i][1] = T
	return new_status,history

def stat(status):
	infe = 0
	for i in status:
		if status[i] != 'S':
			infe += 1
	return infe
	
def statA(status):
	infe = 0
	for i in status:
		if status[i][0] == 'E' or status[i][0] == 'I' or status[i][0] == 'A':
			if status[i] != 'IR' and status[i] != 'AR':
				infe += 1
	return infe

def get_out(name,G,status):
	detect = []
	one_hot = np.identity(len(G))
	Neighbor = set()
	fileobj = open(r'./spreading_data/train', 'w')
	for i in status:
		if status[i][0] == 'I' or status[i] == 'IR':
			fileobj.write(str(i)+'\t')
			for j in one_hot[i-1]:
				fileobj.write(str(j)+'\t')
			fileobj.write(str(1)+'\n')
	fileobj.close()	
