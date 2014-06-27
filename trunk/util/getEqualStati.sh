if [ $# -lt 2 ]
then
    echo 'Usage: tagfile_sortedLen TFfile'
    exit
fi

num=`wc -l $1|cut -f1 -d' '`
one=`echo $num/10|bc`

thre=(0 1 2 3)
suffix=(a b c d e f g h i j)
#suffix=(j i h g f e d c b a)

split -l$one $1 splited

tmpRead=temp$RANDOM.bed
tmpSuffix=tmp$RANDOM
tmpPk=pilePeak$tmpSuffix.bed

for i in ${suffix[*]}
do
    echo $i
    nowF=spliteda$i
    a=`head -n1 $nowF|cut -f5`
    b=`tail -n1 $nowF|cut -f5`
    echo $a $b
    sortBed -i $nowF >$tmpRead
    for j in ${thre[*]}
    do
    	sh callPileUp.sh $tmpRead $j $tmpSuffix
	m=`wc -l $tmpPk|cut -d' ' -f1`
	n=`windowBed -b $tmpPk -a $2 -u -w 0|wc -l`
#	n=`windowBed -a pilePeaktmp.bed -b ../allGm12878TFrmdup.bed -u -w 0|wc -l`
#	echo $n/$m'('`echo $n/$m|bc -l`')'
	echo $n
    done
done
rm $tmpRead
rm $tmpPk
rm spliteda*
