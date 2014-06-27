if [ $# -lt 2 ]
then
    echo 'Usage: bedfile cell'
    exit
fi

echo $1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeq$2CytosolPapMinusRawRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeq$2CytosolPapPlusRawRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeq$2CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeq$2CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqK562R2x75Th1014Il200SigRep1V4.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
