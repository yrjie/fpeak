if [ $# -lt 1 ]
then
    echo 'Usage: sh runTest.sh prefix'
    exit
fi
# $1 prefix
#python callImba.py data/Gm12878Sp1By79_c4.sig data/Gm12878Sp1By79_c4.phyl data/Gm12878Sp1By79_c3.sig data/Gm12878Sp1By79_c3.phyl $1.sig $1.phyl
python callImba.py data/Gm12878FaireSp1_c2.sig data/Gm12878FaireSp1_c2.phyl data/Gm12878FaireSp1_c1.sig data/Gm12878FaireSp1_c1.phyl $1.sig $1.phyl
echo 'count'
#bedfile='data/Gm12878Sp1.bed'
#bedfile='data/Gm12878Sp1By79_c4.bed'
#bedfile='data/Gm12878Sp1By79_c3.bed'
#bedfile='data/Gm12878Sp1_noTr.bed'
bedfile=$1.bed
sort test.lst |uniq -c
paste $bedfile test.lst |awk '{if ($NF==0) print $0}'|cut -f1-4 >temp0
paste $bedfile test.lst |awk '{if ($NF==1) print $0}'|cut -f1-4 >temp1
paste $bedfile test.lst |awk '{if ($NF==2) print $0}'|cut -f1-4 >temp2
echo 'temp0'
closestBed -a temp0 -b ../data/knownGene.bed -t first |cut -f10|sort |uniq -c
echo 'Pol2'
windowBed -a temp0 -b ~/DNase/data/wgEncodeAwgTfbsHaibGm12878Pol24h8Pcr1xUniPk.narrowPeak -u -w 0|wc -l
echo 'RNA-seq'
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig temp0 2 |awk 'BEGIN {l=0;r=0}{if ($1>2*$2) l++; else if ($2>2*$1) r++}END {print l,r}'
sh checkRNA.sh temp0
echo 'temp1'
closestBed -a temp1 -b ../data/knownGene.bed -t first |cut -f10|sort |uniq -c
windowBed -a temp1 -b ~/DNase/data/wgEncodeAwgTfbsHaibGm12878Pol24h8Pcr1xUniPk.narrowPeak -u -w 0|wc -l
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig temp1 2 |awk 'BEGIN {l=0;r=0}{if ($1>2*$2) l++;else if ($2>2*$1) r++}END {print l,r}'
sh checkRNA.sh temp1
echo 'temp2'
closestBed -a temp2 -b ../data/knownGene.bed -t first |cut -f10|sort |uniq -c
windowBed -a temp2 -b ~/DNase/data/wgEncodeAwgTfbsHaibGm12878Pol24h8Pcr1xUniPk.narrowPeak -u -w 0|wc -l
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig temp2 2 |awk 'BEGIN {l=0;r=0}{if ($1>2*$2) l++;else if ($2>2*$1) r++}END {print l,r}'
sh checkRNA.sh temp2
