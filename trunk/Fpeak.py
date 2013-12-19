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
bgWin=20000
plusSig=prefix+'Plus.sig'
minusSig=prefix+'Minus.sig'
plusBed=prefix+'plus.bed'
minusBed=prefix+'minus.bed'
tempBed='temp.bed'
allBin=[]
def genSigFile(infile, outfile):
    fi=open(infile)
    nowChr='chr0'
    nowS=1
    nowE=100
    sum=0
    bp=0
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
	    	#fo.write(('%s\t%d\t%d\t%f\n')%(nowChr,nowS,nowS+bin,1.0*sum/bin))
	    nowChr=chr
	    nowS=start
	    sum=0
	    bp=0
	nowE=end
	sum+=num*(end-start)
	bp+=end-start
    fi.close()
    fo=open(outfile,'w')
    left=right=0
    sum=0
    num=len(allBin)
    print infile+"\t"+str(num)
    for i in range(num):
	while left<i:
	    if allBin[i].beg-allBin[left].beg<=bgWin:
	    	break
	    sum-=allBin[left].val
	    left+=1
	while right<num-1:
	    if allBin[right].end-allBin[i].end>bgWin or allBin[right].chr!=allBin[i].chr:
	    	break
	    right+=1
	    sum+=allBin[right].val
    	allBin[i].val/sum
	fo.write(('%s\t%d\t%d\t%f\n')%(allBin[i].chr,allBin[i].beg,allBin[i].end,2.0*allBin[i].val*bgWin/bin/sum))
    fo.close()

genSigFile(sys.argv[1],plusSig)
genSigFile(sys.argv[2],minusSig)
