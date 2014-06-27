if [ $# -lt 2 ]
then
    echo 'Usage: peakfile uncovTF'
    exit
fi

shuffleBed -i $2 -g /home/ruijie/genome/human.hg19.genome >../shuffledUncovTF.bed
wc -l $1
windowBed -b $2 -a $1 -u -w 0|wc -l
windowBed -a $2 -b $1 -u -w 0|wc -l
windowBed -a ../shuffledUncovTF.bed -b $1 -u -w 0|wc -l
