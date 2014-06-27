import sys,os

if len(sys.argv)<3:
    print 'Usage: bedfile genome'
    exit(1)

weight=[(104,1),
	(164,2),
	(1000,1)]

chrom={}

def getChrLen():
    fi=open(sys.argv[2])
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split()
	chrom[temp[0]]=int(temp[1])
    fi.close()

def readBed():
    nowChr=''
    fi=open(sys.argv[1])
    clen=0
    beg=[0]
    end=[0]
    chrTail=0
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	if temp[0]!=nowChr:
	    if nowChr!='':
	    	printCov(nowChr, chrTail, beg, end)
	    nowChr=temp[0]
	    chrTail=0
	    clen=chrom[nowChr]
	    beg=[0]*(clen+1)
	    end=[0]*(clen+1)
	begF=int(temp[1])
	endF=int(temp[2])
	if endF>chrTail:
	    chrTail=endF
	if begF>clen or endF>clen:
	    raise Exception('Position outside the chromosome: '+line)
	flen=int(temp[4])
	if flen<0:
	    flen=-flen
	w=0
	for x in weight:
	    if flen<x[0]:
	    	w=x[1]
	    	break
	w=1
	if w==0:
	    continue
	beg[begF]+=w
	end[endF]+=w
    fi.close()
    printCov(nowChr, chrTail, beg, end)

def printCov(chr, chrTail, beg, end):
    pre=0
    sig=0
    for i in range(chrTail+1):
    	if beg[i]==0 and end[i]==0:
    	    continue
	if beg[i]!=end[i]:
	    if sig>0:
	    	print '\t'.join([chr, str(pre), str(i), str(sig)])
	    pre=i
	sig+=beg[i]
	sig-=end[i]

if __name__ =='__main__':
    getChrLen()
    readBed()
