import sys,os
import scipy.stats as ss
import scipy as sp
import time
import bisect
import math

if len(sys.argv)<2:
    print 'Usage: python fitGamma.py infile'
    exit(1)

lines=open(sys.argv[1]).readlines()
data=[float(x) for x in lines ]#if float(x)>0]
#data=ss.gamma.rvs(alpha,loc=loc,scale=beta,size=10000)

a=ss.gamma.fit(data)
fit_alpha,fit_loc,fit_beta=ss.gamma.fit(data)
#(fit_loc,fit_scale)=ss.cauchy.fit(data)
#print sys.argv[1]
#print a
#print(fit_loc,fit_scale)
#print ss.cauchy.pdf(2.9, loc=fit_loc,scale=fit_scale)
#print ss.gamma.pdf(2.9, a[0],loc=a[1],scale=a[2])
beg=time.time()
[left,right]=ss.gamma.interval(0.999,a[0],loc=a[1],scale=a[2])
#print left,right
x=[left+i*(right-left)/1000.0 for i in range(1001)]
y=[ss.gamma.pdf(i, a[0],loc=a[1],scale=a[2]) for i in x]
#z=[sys.argv[1]]+[str(i) for i in x]
z=ss.gamma.rvs(a[0],loc=a[1],scale=a[2],size=10000)
print "\n".join([str(i) for i in z])
ind=[bisect.bisect_left(x,i,lo=0,hi=len(x)-1) for i in data]
res1=[y[i] for i in ind]
#res2=[ss.gamma.pdf(i, a[0],loc=a[1],scale=a[2]) for i in data]
res2=[ss.gamma.pdf(i, a[0],loc=a[1],scale=a[2]) for i in x]
#print sum(res2)*(right-left)
#res=sp.mean([abs(math.log(res1[i])-math.log(res2[i])) for i in range(len(res1))])
#res=sp.mean([abs(res1[i]-res2[i]) for i in range(len(res1))])
#print res
#print int((time.time()-beg)*1000)
#print "cdf: "+str(ss.gamma.cdf(10, a[0],loc=a[1],scale=a[2]))
#print ss.norm.fit(data)
#print ss.gamma.sf(10,1.763804808584083e-05, loc=2.9411799999999996e-09, scale=0.34638959197470232)

#print ss.gamma.sf(50,fit_alpha,loc=fit_loc,scale=fit_beta)
#print ss.gamma.interval(0.9,fit_alpha,loc=fit_loc,scale=fit_beta)
#print ss.cauchy.interval(0.9,loc=fit_loc,scale=fit_scale)
