if [ $# -lt 2 ]
then
    echo 'Usage: pos.tab neg.tab'
    exit
fi

posF=$1
negF=$2

rm temp.dat
for ((i=1;i<=5;i++)) do
    cut -f $i $posF|awk -v name=`echo $i*100|bc` '{if (NR==1) printf name; printf "\t"$1}END{printf "\n"}' >>temp.dat
done
R --no-save --slave --args temp.dat <~/bin/plotBox.R
mv temp.png pos.png

rm temp.dat
for ((i=1;i<=5;i++)) do
    cut -f $i $negF|awk -v name=`echo $i*100|bc` '{if (NR==1) printf name; printf "\t"$1}END{printf "\n"}' >>temp.dat
done
R --no-save --slave --args temp.dat <~/bin/plotBox.R
mv temp.png neg.png

