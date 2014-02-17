import sys,os
import scipy.stats as ss
import scipy as sp

if len(sys.argv)<2:
    print 'Usage: python fitGamma.py infile'
    exit(1)

lines=open(sys.argv[1]).readlines()
data=[float(x) for x in lines if float(x)>0]
#data=ss.gamma.rvs(alpha,loc=loc,scale=beta,size=10000)

fit_alpha,fit_loc,fit_beta=ss.gamma.fit(data)
print sys.argv[1]
print(fit_alpha,fit_loc,fit_beta)
#print ss.gamma.sf(10,1.763804808584083e-05, loc=2.9411799999999996e-09, scale=0.34638959197470232)

#print ss.gamma.sf(50,fit_alpha,loc=fit_loc,scale=fit_beta)
print ss.gamma.interval(0.9,fit_alpha,loc=fit_loc,scale=fit_beta)
