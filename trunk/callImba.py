import os,sys
import numpy as np
import time
from hmmpytk import hmm_faster

if len(sys.argv)<3:
    print 'Usage: python callImba.py FaireSig conservation'
    exit(1)

allObs=[]
meanObs=[]

def readData():
    rowSig=[]
    rowCons=[]
    sigFile=open(sys.argv[1])
    consFile=open(sys.argv[2])
    for sig in sigFile:
	sig=sig.strip()
	cons=consFile.readline().strip()
	if len(sig)<1 or len(cons)<1:
	    continue
	tempSig=[float(x) for x in sig.split('\t')]
	meanSig=np.mean(tempSig)
	stdSig=np.std(tempSig)
	#tempCons=[float(x) if float(x)>0 else 0 for x in cons.split('\t')]
	tempCons=[float(x) for x in cons.split('\t')]
	meanCons=np.mean(tempCons)
	stdCons=np.std(tempCons)
	obs=[]
	if len(rowSig)<1:
	    rowSig=[0]*len(tempSig)
	    rowCons=[0]*len(tempCons)
	for i in xrange(len(tempSig)):
	    state=''
	    if tempSig[i]-meanSig>stdSig:
		state+='hf'
	    elif tempSig[i]-meanSig<-0.5*stdSig:
		state+='lf'
	    else:
		state+='mf'
	    if tempCons[i]-meanCons>stdCons:
		state+='hc'
	    elif tempCons[i]-meanCons<-0.5*stdCons:
		state+='lc'
	    else:
		state+='mc'
	    obs.append(state)
	    rowSig[i]+=tempSig[i]
	    rowCons[i]+=tempCons[i]
	allObs.append(obs)
    sigFile.close()
    consFile.close()
    meanSig=np.mean(rowSig)
    stdSig=np.std(rowSig)
    meanCons=np.mean(rowCons)
    stdCons=np.std(rowCons)
    for i in xrange(len(rowSig)):
    	state=''
	if rowSig[i]-meanSig>stdSig:
	    state+='hf'
	elif rowSig[i]-meanSig<-0.5*stdSig:
	    state+='lf'
	else:
	    state+='mf'
	if rowCons[i]-meanCons>stdCons:
	    state+='hc'
	elif rowCons[i]-meanCons<-0.5*stdCons:
	    state+='lc'
	else:
	    state+='mc'
	meanObs.append(state)

def initMat(Pi, E, T):
    for x in Pi:
    	Pi[x]=0
    for x in E:
	for y in E[x]:
	    E[x][y]=0
    for x in T:
	for y in T[x]:
	    T[x][y]=0

def runHmm():
    hmm_model = hmm_faster.HMM()
    #hmm_model.set_states(['Bg','Peak','Active', 'Drop'])
    #hmm_model.set_states(['Peak','Active', 'Drop','Bg'])
    hmm_model.set_states(['Active', 'Peak', 'Drop','Bg'])
    hmm_model.set_observations(['hfhc', 'hfmc', 'hflc', 'mfhc', 'mfmc', 'mflc','lfhc', 'lfmc', 'lflc'])
    hmm_model.randomize_matrices(seed = int(time.time()))
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 7, 'Bg':3}
    T_matrix={'Peak':{'Peak': 0.3, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.4},
	    'Drop':{'Peak': 0.5, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.1},
	    'Active':{'Peak': 0.1, 'Drop': 0.3, 'Active': 0.3, 'Bg':0.1},
	    'Bg':{'Peak': 0.2, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.7}}
    E_matrix={'Peak':{'hfhc': 0.4, 'hfmc': 0.4, 'hflc': 0.2,'mfhc': 0, 'mfmc': 0, 'mflc': 0,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
	    'Drop':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.1, 'mfmc': 0.1, 'mflc': 0,'lfhc': 0.5, 'lfmc': 0.48, 'lflc': 0.02},
	    'Active':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.3, 'mfmc': 0.5, 'mflc': 0.2,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
	    'Bg':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0, 'mfmc': 2, 'mflc': 3,'lfhc': 0, 'lfmc': 2, 'lflc': 3}}
    for i in range(10):
	hmm_model.set_initial_matrix(Pi_matrix)
	hmm_model.set_transition_matrix(T_matrix)
	hmm_model.set_emission_matrix(E_matrix)
	#print hmm_model.get_initial_matrix()
	hmm_model.train(meanObs, max_iteration=1000, delta=0.001)
	#print hmm_model.get_model()
	left=0
	right=0
	mid=0
	#print meanObs
	#print hmm_model.viterbi(meanObs)
	initMat(Pi_matrix, E_matrix, T_matrix)
	for x in allObs:
	    result=hmm_model.viterbi(x)
	    #print x
	    #print result
	    Pi_matrix[result[0]]+=1
	    for st in Pi_matrix:
	    	lstSt=[i for i in range(len(result)) if result[i]==st]
		for j in lstSt:
		    y=x[j]
		    E_matrix[st][y]+=1
		for j in lstSt:
		    if j+1>=len(result):
		    	continue
		    y=result[j+1]
		    T_matrix[st][y]+=1
	    if 'Peak' not in result:
		mid+=1
		continue
	    peakPos=int(np.median([i for i in xrange(len(result)) if result[i]=='Peak']))
	    lcnt=len([i for i in xrange(0,peakPos) if result[i]=='Drop'])
	    rcnt=len([i for i in xrange(peakPos,len(result)) if result[i]=='Drop'])
	    #lcnt=len([i for i in xrange(0,peakPos) if result[i]=='Drop' or result[i]=='Active'])
	    #rcnt=len([i for i in xrange(peakPos,len(result)) if result[i]=='Drop' or result[i]=='Active'])
	    if lcnt-rcnt>0:
		left+=1
	    elif lcnt-rcnt<0:
		right+=1
	    else:
		mid+=1
	print str(left)+"\t"+str(mid)+"\t"+str(right)
    #hmm_model.train(['high','mid','high','low','high','mid','high','low','high'], max_iteration=100, delta=0.0001)
    #print hmm_model.viterbi(['low', 'high', 'high', 'low'])
    #print hmm_model.get_model()

readData()
runHmm()
