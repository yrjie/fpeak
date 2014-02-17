import os,sys
import numpy as np
from hmmpytk import hmm_faster

if len(sys.argv)<3:
    print 'Usage: python callImba.py FaireSig conservation'
    exit(1)

allSig=[]
sigFile=open(sys.argv[1])
for line in sigFile:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(x) for x in line.split('\t')]
    meanX=np.mean(temp)
    obs=[]
    for x in temp:
	if x/meanX>2:
	    obs.append('high')
	elif x/meanX<0.5:
	    obs.append('low')
	else:
	    obs.append('mid')
    allSig.append(obs)

allCons=[]
consFile=open(sys.argv[2])
for line in consFile:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(x) if float(x)>0 else 0 for x in line.split('\t')]
    meanX=np.mean(temp)
    obs=[]
    for x in temp:
	if x/meanX>2:
	    obs.append('high')
	elif x/meanX<0.5:
	    obs.append('low')
	else:
	    obs.append('mid')
    allCons.append(obs)

hmm_model = hmm_faster.HMM()
hmm_model.set_states(['Peak', 'Drop', 'Active', 'Bg'])
hmm_model.set_observations(['high', 'mid', 'low'])
hmm_model.randomize_matrices(seed = 100)
Pi_matrix={'Peak': 0, 'Drop': 0, 'Active': 0, 'Bg':1}
T_matrix={'Peak':{'Peak': 0.4, 'Drop': 0.3, 'Active': 0.2, 'Bg':0.1},
	'Drop':{'Peak': 0.4, 'Drop': 0.3, 'Active': 0.2, 'Bg':0.1},
	'Active':{'Peak': 0.3, 'Drop': 0.3, 'Active': 0.2, 'Bg':0.1},
	'Bg':{'Peak': 0.2, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.5}}
E_matrix={'Peak':{'high': 0.9, 'mid': 0.08, 'low': 0.02},
	'Drop':{'high': 0.2, 'mid': 0.3, 'low': 0.5},
	'Active':{'high': 0.2, 'mid': 0.5, 'low': 0.3},
	'Bg':{'high': 0.1, 'mid': 0.8, 'low': 0.1}}
hmm_model.set_initial_matrix(Pi_matrix)
#hmm_model.set_transition_matrix(T_matrix)
#hmm_model.set_emission_matrix(E_matrix)
hmm_model.train(['high','mid','high','low','high','mid','high','low','high'], max_iteration=100, delta=0.0001)
for x in allSig:
    print hmm_model.viterbi(x)
#print hmm_model.get_model()
#hmm_model.train(['high','mid','high','low','high','mid','high','low','high'], max_iteration=100, delta=0.0001)
#print hmm_model.viterbi(['low', 'high', 'high', 'low'])
#print hmm_model.get_model()
