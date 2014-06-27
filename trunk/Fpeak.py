import sys,os

if len(sys.argv)<4:
    print 'Usage: python Fpeak.py plus.bedGraph minus.bedGraph prefix'
    exit(1)

class ChrBin:
    chr='chr1'
    beg=0
    end=0
    val=0
    def __init__(self, chr, beg, end, val):
    	self.chr=chr
    	self.beg=beg
    	self.end=end
    	self.val=val

prefix=sys.argv[3]
bin=200
bgWin=2000
plusSig=prefix+'Plus.sig'
minusSig=prefix+'Minus.sig'
plusBed=prefix+'plus.bed'
minusBed=prefix+'minus.bed'
tempBed='temp.bed'

def genSigFile(infile, outfile):
    allBin=[]
    allSum=[]
    fi=open(infile)
    fo=open(outfile,'w')
    nowChr='chr0'
    nowS=1
    nowE=100
    smtP=0
    smtSig=0
    sum=0
    bp=0
    csum=0
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	chr=temp[0]
	start=int(temp[1])
	end=int(temp[2])
	num=int(temp[3])
	if chr!=nowChr or end-nowS>bin:
	    if nowChr!='chr0':
	    	allBin.append(ChrBin(nowChr, nowS, nowS+bin, sum))
	    	csum+=sum
	    	allSum.append(csum)
		if 'Plus' in outfile:
		    nowE=nowS+bin
		else:
		    nowS=max(1,nowE-bin)
	    	fo.write(('%s\t%d\t%d\t%d\t%f\n')%(nowChr,nowS,nowE,smtP,1.0*sum/bin))
	    nowChr=chr
	    nowS=start
	    sum=0
	    smtP=0
	    smtSig=0
	    bp=0
	nowE=end
	sum+=num*(end-start)
	if num>smtSig:
	    smtP=(start+end)/2
	    smtSig=num
	bp+=end-start
    fi.close()
    fo.close()
#    fo=open(outfile,'w')
#    left=right=0
#    num=len(allBin)
#    print infile+"\t"+str(num)
#    for i in range(num):
#	while left<i:
#	    if allBin[i].chr==allBin[left].chr and allBin[i].beg-allBin[left].beg<=bgWin:
#	    	break
#	    left+=1
#	while right<num-1:
#	    if allBin[right].end-allBin[i].end>bgWin or allBin[right].chr!=allBin[i].chr:
#	    	break
#	    right+=1
#	if left==0:
#	    sum=allSum[right-1]
#	else:
#	    sum=allSum[right-1]-allSum[left]
#	if sum==0:
#	    sum=1
#	fo.write(('%s\t%d\t%d\t%f\n')%(allBin[i].chr,allBin[i].beg,allBin[i].end,2.0*allBin[i].val*bgWin/bin/sum))
#	#fo.write(('%s\t%d\t%d\t%f\n')%(allBin[i].chr,allBin[i].beg,allBin[i].end,1.0*allBin[i].val/bin))
#    fo.close()

def refinePeak(peakfile, plustab, minustab):
    fi=open(peakfile)
    fiP=open(plustab)
    fiM=open(minustab)
    fo=open(prefix+'Peak.bed','w')
    for peak in fi:
    	peak=peak.strip()
	if len(peak)<1:
	    continue

genSigFile(sys.argv[1],plusSig)
genSigFile(sys.argv[2],minusSig)
#cmd='sh mergePeak.sh '+plusSig+' '+minusSig+' '+prefix
#os.system(cmd)
#refinePeak(prefix+'.bed.temp','tempPlus.tab','tempMinus.tab')
