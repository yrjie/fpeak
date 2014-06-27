if [ $# -lt 2 ]
then
    echo 'Usage: peakfile oneHead'
    exit
fi

#one=20000
one=$2

for ((i=1;i<=10;i++)) do
    num=`echo $i*$one|bc`
    sort -k4gr $1 |head -n $num|windowBed -a stdin -b ../ATAC/uncovTFrmdup_pile.bed -u -w 0|wc -l
    #sort -k4gr $1 |head -n $num|windowBed -b stdin -a ../ATAC/uncovTFrmdup_pile.bed -u -w 0|wc -l
done

