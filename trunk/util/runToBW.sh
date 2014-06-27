if [ $# -lt 1 ]
then
    echo "Usage: bamfile"
    exit
fi
# $1 input bam
bedF=${1%\.*}.bed
python ~/bin/extendTag.py $1 ~/genome/mouse.mm9.genome 150 >$bedF
num=`wc -l $bedF |cut -d' ' -f1`
echo $bedF $num
sc=`echo 1000000/$num |bc -l`
file=${1##*/}
bedtools genomecov -i $bedF -g ~/genome/mouse.mm9.genome -bg -scale $sc >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/mouse.mm9.genome coverage/${file%\.*}.bigWig
#bigWigSummaryBatch 2bp/$head.minus.bigWig data/wgEncodeAwgTfbsUtaK562CtcfUniPk_plus.narrowPeak 1 -type=max >pattern/K562Ctcf_plus.xls
