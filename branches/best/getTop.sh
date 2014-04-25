sort -gr -k4 H1hescPeak.bed |head -n400000|cut -f1-3|centerBed |slopBed -g ~/genome/human.hg19.genome -b 500 >temp
paste temp tempall.sig.prob |sort -k5 -gr|head -n100000 >temp1
