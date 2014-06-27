if [ $# -lt 4 ]
then
    echo "Usage: sh getRocPeak.sh peakfile cell TFfile suffix"
    exit
fi

file=$1
cell=$2
TFfile=$3
suffix=Peak$4
rm output/$cell$suffix.tab
num=`wc -l $file |cut -d' ' -f1`
one=`echo $num/10 |bc`

echo $cell$suffix.tab
#max1=10000 # for ATAC under DFilter
#max1=20000 # for ATAC
max1=50000 # for ATAC_low
#max1=100000 # for Faire
#max1=200000 # for Faire

#if [ $max1 -lt $one ]
#then
#    one=$max1
#fi
sortFile=sorted$RANDOM
# be careful about the column of score
sort -k4 -gr $file >$sortFile
numP=`windowBed -a $sortFile -b $TFfile -u -w 0 |wc -l`
numN=`windowBed -a $sortFile -b $TFfile -v -w 0 |wc -l`
for ((i=1;i<=10;i++)) do
    now=`echo $one*$i|bc`
    head -n $now $sortFile >tempRoc$suffix
    nowP=`windowBed -a tempRoc$suffix -b $TFfile -u -w 0 |wc -l`
    nowN=`windowBed -a tempRoc$suffix -b $TFfile -v -w 0 |wc -l`
    tpr=`echo $nowP/$numP|bc -l`
    fpr=`echo $nowN/$numN|bc -l`
    echo -e $fpr"\t"$tpr >>output/$cell$suffix.tab
#    sh checkFpCov.sh $cell tempRoc$suffix $suffix
#    awk 'BEGIN {FS=OFS="\t";fp=0;cov=0} {fp+=$5;cov+=$4} END {print fp/NR,cov/NR}' output/$cell$suffix'_all'.xls >>output/$cell$suffix.tab
done
rm $sortFile
