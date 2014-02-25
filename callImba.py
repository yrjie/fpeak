import os,sys
import numpy as np
import time
import scipy.stats as ss
from hmmpytk import hmm_faster
from hmmpytk import hmm_gamma

if len(sys.argv)<7:
    print 'Usage: python callImba.py FaireSig conservation testFaire testCons'
    exit(1)

# class RealBinObs:
#     fsig=0
#     cons=0
#     def __init__(self, f, c):
#         self.fsig=f
#         self.cons=c

realL=[]
realR=[]
realAll=[]
trainLObs=[]
trainRObs=[]
allObs=[]
meanObs=[]

def readData(obsLst, realLst, fileF, fileC):
    rowSig=[]
    rowCons=[]
    sigFile=open(fileF)
    consFile=open(fileC)
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
        real=[]
        if len(rowSig)<1:
            rowSig=[0]*len(tempSig)
            rowCons=[0]*len(tempCons)
        for i in xrange(len(tempSig)):
            real.append([tempSig[i], tempCons[i]])
#             try:
#                 real.append(RealBinObs((tempSig[i]-meanSig)/stdSig, (tempCons[i]-meanCons)/stdCons))
#             except:
#                 real.append(RealBinObs(0, 0))
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
        obsLst.append(obs)
        realLst.append(real)
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

def buildLeft():
    hmm_model = hmm_faster.HMM()
    hmm_model.set_states(['Active', 'Peak', 'Drop','Bg'])
    hmm_model.set_observations(['hfhc', 'hfmc', 'hflc', 'mfhc', 'mfmc', 'mflc','lfhc', 'lfmc', 'lflc'])
    hmm_model.randomize_matrices(seed = int(time.time()))
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 7, 'Bg':3}
    T_matrix={'Peak':{'Peak': 0.3, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.4},
        'Drop':{'Peak': 0.5, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.1},
        'Active':{'Peak': 0.1, 'Drop': 0.3, 'Active': 0.3, 'Bg':0.1},
        'Bg':{'Peak': 0.2, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.7}}
    E_matrix={'Peak':{'hfhc': 0.4, 'hfmc': 0.4, 'hflc': 0.2,'mfhc': 0.4, 'mfmc': 0.4, 'mflc': 0.1,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Drop':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.1, 'mfmc': 0.1, 'mflc': 0,'lfhc': 0.5, 'lfmc': 0.48, 'lflc': 0.02},
        'Active':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.3, 'mfmc': 0.5, 'mflc': 0.2,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Bg':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0, 'mfmc': 2, 'mflc': 3,'lfhc': 0, 'lfmc': 2, 'lflc': 3}}
    hmm_model.set_initial_matrix(Pi_matrix)
    hmm_model.set_transition_matrix(T_matrix)
    hmm_model.set_emission_matrix(E_matrix)
    #print hmm_model.get_initial_matrix()
    hmm_model.train(sum(trainLObs,[]), max_iteration=1000, delta=0.001)
    return hmm_model

def buildMid():
    hmm_model = hmm_faster.HMM()
    hmm_model.set_states(['Active', 'Peak', 'Drop','Bg'])
    hmm_model.set_observations(['hfhc', 'hfmc', 'hflc', 'mfhc', 'mfmc', 'mflc','lfhc', 'lfmc', 'lflc'])
    hmm_model.randomize_matrices(seed = int(time.time()))
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 5, 'Bg':5}
    T_matrix={'Peak':{'Peak': 0.3, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.4},
        'Drop':{'Peak': 0.5, 'Drop': 0.1, 'Active': 0.5, 'Bg':0.1},
        'Active':{'Peak': 0.1, 'Drop': 0.3, 'Active': 0.3, 'Bg':0.1},
        'Bg':{'Peak': 0.4, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.4}}
    E_matrix={'Peak':{'hfhc': 0.4, 'hfmc': 0.4, 'hflc': 0.2,'mfhc': 0, 'mfmc': 0, 'mflc': 0,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Drop':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.1, 'mfmc': 0.1, 'mflc': 0,'lfhc': 0.5, 'lfmc': 0.48, 'lflc': 0.02},
        'Active':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.3, 'mfmc': 0.5, 'mflc': 0.2,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Bg':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0, 'mfmc': 2, 'mflc': 3,'lfhc': 0, 'lfmc': 2, 'lflc': 3}}
    hmm_model.set_initial_matrix(Pi_matrix)
    hmm_model.set_transition_matrix(T_matrix)
    hmm_model.set_emission_matrix(E_matrix)
    #print hmm_model.get_initial_matrix()
    return hmm_model

def buildRight():
    hmm_model = hmm_faster.HMM()
    hmm_model.set_states(['Active', 'Bg', 'Drop', 'Peak'])
    hmm_model.set_observations(['hfhc', 'hfmc', 'hflc', 'mfhc', 'mfmc', 'mflc','lfhc', 'lfmc', 'lflc'])
    hmm_model.randomize_matrices(seed = int(time.time()))
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 3, 'Bg':7}
    T_matrix={'Peak':{'Peak': 0.6, 'Drop': 0.5, 'Active': 0.1, 'Bg':0.1},
        'Drop':{'Peak': 0.1, 'Drop': 0.2, 'Active': 0.5, 'Bg':0.1},
        'Active':{'Peak': 0.1, 'Drop': 0.1, 'Active': 0.6, 'Bg':0.1},
        'Bg':{'Peak': 0.5, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.7}}
    E_matrix={'Peak':{'hfhc': 0.4, 'hfmc': 0.4, 'hflc': 0.2,'mfhc': 0, 'mfmc': 0, 'mflc': 0,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Drop':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.1, 'mfmc': 0.1, 'mflc': 0,'lfhc': 0.5, 'lfmc': 0.48, 'lflc': 0.02},
        'Active':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0.3, 'mfmc': 0.5, 'mflc': 0.2,'lfhc': 0, 'lfmc': 0, 'lflc': 0},
        'Bg':{'hfhc': 0, 'hfmc': 0, 'hflc': 0,'mfhc': 0, 'mfmc': 2, 'mflc': 3,'lfhc': 0, 'lfmc': 2, 'lflc': 3}}
    hmm_model.set_initial_matrix(Pi_matrix)
    hmm_model.set_transition_matrix(T_matrix)
    hmm_model.set_emission_matrix(E_matrix)
    #print hmm_model.get_initial_matrix()
    hmm_model.train(sum(trainRObs, []), max_iteration=1000, delta=0.001)
    return hmm_model

def buildLeftGamma(gamma_par):
    hmm_model = hmm_gamma.HMMgamma()
    hmm_model.set_states(['Active', 'Peak', 'Drop','Bg'])
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 7, 'Bg':3}
    T_matrix={'Peak':{'Peak': 0.3, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.4},
        'Drop':{'Peak': 0.5, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.1},
        'Active':{'Peak': 0.1, 'Drop': 0.3, 'Active': 0.3, 'Bg':0.1},
        'Bg':{'Peak': 0.2, 'Drop': 0.2, 'Active': 0.1, 'Bg':0.7}}
    hmm_model.set_initial_matrix(Pi_matrix)
    hmm_model.set_transition_matrix(T_matrix)
    hmm_model.set_emission_table(gamma_par)
#     hmm_model.train(sum(realL,[]), max_iteration=10, delta=0.001)
    return hmm_model

def buildRightGamma(gamma_par):
    hmm_model = hmm_gamma.HMMgamma()
    hmm_model.set_states(['Active', 'Peak', 'Drop','Bg'])
    Pi_matrix={'Peak': 1, 'Drop': 1, 'Active': 7, 'Bg':3}
    T_matrix={'Peak':{'Peak': 0.6, 'Drop': 0.5, 'Active': 0.1, 'Bg':0.1},
        'Drop':{'Peak': 0.1, 'Drop': 0.2, 'Active': 0.5, 'Bg':0.1},
        'Active':{'Peak': 0.1, 'Drop': 0.1, 'Active': 0.6, 'Bg':0.1},
        'Bg':{'Peak': 0.5, 'Drop': 0.1, 'Active': 0.1, 'Bg':0.7}}
    hmm_model.set_initial_matrix(Pi_matrix)
    hmm_model.set_transition_matrix(T_matrix)
    hmm_model.set_emission_table(gamma_par)
    #print hmm_model.get_initial_matrix()
#     hmm_model.train(sum(realR,[]), max_iteration=10, delta=0.001)
    return hmm_model

def getDistrByState(realLst, obsLst, hmm_model):
    states=['Active','Bg','Drop','Peak']
    realSig=[]
    realCons=[]
    gamma_par={}
#     result=hmm_model.viterbi(meanObs)
#     print meanObs
#     print result
    for st in states:
        realSig.append([])
        realCons.append([])
    for i in xrange(len(obsLst)):
        x=obsLst[i]
        result=hmm_model.viterbi(x)
        for j in xrange(len(result)):
            re=result[j]
            for k in xrange(len(states)):
                if re==states[k]:
                    realSig[k].append(realLst[i][j][0]+0.1)
                    realCons[k].append(realLst[i][j][1])
                    break
    for i in xrange(len(states)):
        par_for_st=[]
        par_for_st.append(ss.gamma.fit(realSig[i]))
        par_for_st.append(ss.gamma.fit(realCons[i]))
        gamma_par[states[i]]=par_for_st
    return gamma_par
    
def runHmm():
    hmm_all=[]
    beg=time.time()
    hmm_all.append(buildLeft())
    print 1
    gamma_parL=getDistrByState(realL, trainLObs, hmm_all[0])
    print 2
#     hmm_all[0]=buildLeftGamma(gamma_parL)
    print 3
    hmm_all.append(buildMid())
    print 4
    hmm_all.append(buildRight())
    print 5
    gamma_parR=getDistrByState(realR, trainRObs, hmm_all[2])
    print 6
#     hmm_all[2]=buildRightGamma(gamma_parR)
#     gamma_parR=getDistrByState(realR, trainRObs, hmm_all[2])
    trainTS=time.time()
    print 'trainTS: '+str(trainTS-beg)
    left=0
    right=0
    mid=0
    ind=[0 for i in xrange(len(realAll))]
    for k in xrange(len(realAll)):
#         x=realAll[k]
        x=allObs[k]
        ma=0
        ind[k]=0
# 	result=hmm_model.viterbi(x)
# 	print x
# 	print result
# 	    Pi_matrix[result[0]]+=1
# 	    for st in Pi_matrix:
# 	    	lstSt=[i for i in range(len(result)) if result[i]==st]
# 		for j in lstSt:
# 		    y=x[j]
# 		    E_matrix[st][y]+=1
# 		for j in lstSt:
# 		    if j+1>=len(result):
# 		    	continue
# 		    y=result[j+1]
# 		    T_matrix[st][y]+=1
        pr0=hmm_all[0].evaluate(x)
        pr2=hmm_all[2].evaluate(x)
#         print str(pr0)+"\t"+str(pr2)
        if pr0>2*pr2:
            ind[k]=0
        elif pr2>2*pr0:
            ind[k]=2
        else:
            ind[k]=1
#         for i in xrange(len(hmm_all)):
#             now=hmm_all[i].evaluate(x)
#             if now>ma:
#                 ma=now
#                 ind[k]=i
# 	    peakPos=int(np.median([i for i in xrange(len(result)) if result[i]=='Peak']))
# 	    lcnt=len([i for i in xrange(0,peakPos) if result[i]=='Drop'])
# 	    rcnt=len([i for i in xrange(peakPos,len(result)) if result[i]=='Drop'])
	    #lcnt=len([i for i in xrange(0,peakPos) if result[i]=='Drop' or result[i]=='Active'])
	    #rcnt=len([i for i in xrange(peakPos,len(result)) if result[i]=='Drop' or result[i]=='Active'])
        if ind[k]==0:
            left+=1
        elif ind[k]==2:
            right+=1
        else:
            mid+=1
#     for i in xrange(len(hmm_all)):
#         print hmm_all[i].viterbi(meanObs)
#         print hmm_all[i].evaluate(meanObs)
    testTS=time.time()
    fo=open('test.lst','w')
    fo.write('\n'.join([str(x) for x in ind]))
    fo.close()
#     print str(left)+"\t"+str(mid)+"\t"+str(right)
    #print "Train: "+str(trainTS-beg)+" Test: "+str(testTS-trainTS)
	#print hmm_model.evaluate(meanObs)
    #hmm_model.train(['high','mid','high','low','high','mid','high','low','high'], max_iteration=100, delta=0.0001)
    #print hmm_model.viterbi(['low', 'high', 'high', 'low'])
    #print hmm_model.get_model()

readData(trainLObs, realL, sys.argv[1], sys.argv[2])
readData(trainRObs, realR, sys.argv[3], sys.argv[4])
readData(allObs, realAll, sys.argv[5], sys.argv[6])
runHmm()
