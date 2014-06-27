import sys,os

if len(sys.argv)<2:
    print 'Usage: bedfile'
    exit(1)

sumDict={}
cntDict={}
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[0] not in sumDict:
    	sumDict[temp[0]]=0
    	cntDict[temp[0]]=0
    num=int(temp[2])-int(temp[1])
    sumDict[temp[0]]+=num
    cntDict[temp[0]]+=1
fi.close()
for x in sumDict:
    meanL=1.0*sumDict[x]/cntDict[x]
    print '\t'.join([x,str(meanL)])
