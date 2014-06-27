if [ $# -lt 3 ]
then
    echo "Usage: sig1Peak overPeak windowSize"
    exit
fi

echo $3
a=`windowBed -a $1 -b $2 -u -w $3 |wc -l`
b=`windowBed -a $1 -b $2 -u -w $3 |windowBed -a stdin -b ../allGm12878TFrmdup.bed -u -w 0 |wc -l`
c=`windowBed -a $1 -b $2 -v -w $3 |wc -l`
d=`windowBed -a $1 -b $2 -v -w $3 |windowBed -a stdin -b ../allGm12878TFrmdup.bed -u -w 0 |wc -l`
echo $a $b $c $d
