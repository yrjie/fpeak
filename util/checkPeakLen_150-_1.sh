if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

for ((i=50;i<=300;i+=25)) do
    n1=`awk -v len=$i '{if ($3-$2>len) print $0}' $1 |wc -l`
    n2=`awk -v len=$i '{if ($3-$2>len) print $0}' $1 |windowBed -a stdin -b ../allGm12878TFrmdup.bed -u -w 0|wc -l`
    echo -e $n1"\t"$n2|awk 'BEGIN{OFS="\t"}{print $1,$2,$2/$1}'
done
