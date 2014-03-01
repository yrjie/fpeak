echo $1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapMinusRawRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapPlusRawRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqK562R2x75Th1014Il200SigRep1V4.bigWig $1 1 |awk 'BEGIN {s=0}{s+=$1} END {print s/NR} '
