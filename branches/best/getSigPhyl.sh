if [ $# -lt 2 ]
then
    echo "Usage: sh getSigPhyl.sh prefix FaireFile"
    exit
fi
bigWigSummaryBatch $2 $1.bed 20 >$1.sig
bigWigSummaryBatch ~/data/hg19/phyloP46wayAll.bw $1.bed 20 >$1.phyl
