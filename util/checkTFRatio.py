import sys,os
from subprocess import *
import operator

if len(sys.argv)<3:
    print 'Usage: peakFile TFfile'
    exit(1)

peak=sys.argv[1]
tf=sys.argv[2]

fi=open(tf)
allNum={}
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[3] in allNum:
    	allNum[temp[3]]+=1
    else:
    	allNum[temp[3]]=1
fi.close()

cmd=('windowBed -a '+tf+' -b '+peak+' -u -w 0').split()
p=Popen(cmd,stdout=PIPE)
overlap=p.stdout.readlines()
covTF={}
for line in overlap:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[3] in covTF:
    	covTF[temp[3]]+=1
    else:
    	covTF[temp[3]]=1

for x in covTF:
    covTF[x]/=1.0*allNum[x]

sortedTF = sorted(covTF.iteritems(), key=operator.itemgetter(1),reverse=True)
for x in sortedTF:
    print '\t'.join([str(x[0]),str(x[1])])
