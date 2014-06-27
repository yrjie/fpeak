import sys,os

if len(sys.argv)<2:
    print 'Usage: overlapFile'
    exit(1)

nowID=''
print '\t'.join(['peak','TF','<100','100-200','200-300','300-400','>400','total'])
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    id=temp[0]+':'+temp[1]+'-'+temp[2]+'\t'+temp[3]
    if id!=nowID:
	if nowID!='':
	    print '\t'.join([nowID]+[str(x) for x in distr]+[str(sum(distr))])
    	distr=[0]*5
    	nowID=id
    flen=int(temp[-2])
#    print flen
    if flen<100:
    	distr[0]+=1
    elif flen<200:
    	distr[1]+=1
    elif flen<300:
    	distr[2]+=1
    elif flen<400:
    	distr[3]+=1
    else:
    	distr[4]+=1
fi.close()
