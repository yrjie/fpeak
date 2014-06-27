if [ $# -lt 2 ]
then
    echo 'Usage: tagfile_sorted TFfile'
    exit
fi

len=(0 150 200 250 300 350 400 1000 1000000)
thre=(0 1 2 3)

num=${#len[*]}
tmpRead=temp$RANDOM.bed
tmpSuffix=tmp$RANDOM
tmpPk=pilePeak$tmpSuffix.bed

for ((i=1;i<$num;i++))
do
    echo ${len[i-1]}
    a=${len[i-1]}
    b=${len[i]}
    awk -v a=$a -v b=$b 'BEGIN{OFS="\t"}{if ($5<0) $5=-$5; if ($5>=a&&$5<b) print $0}' $1 >$tmpRead
    for j in ${thre[*]}
    do
	sh callPileUp.sh $tmpRead $j $tmpSuffix
	m=`wc -l $tmpPk|cut -d' ' -f1`
	n=`windowBed -b $tmpPk -a $2 -u -w 0|wc -l`
#	n=`windowBed -a pilePeaktmp.bed -b $2 -u -w 0|wc -l`
#	echo $n/$m'('`echo $n/$m|bc -l`')'
	echo $n
    done
done
rm $tmpRead
rm $tmpPk
