if [ $# -lt 1 ]
then
    echo "bedfileL bedfileM bedfileR"
    exit
fi
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig $1 2 >temp
#paste $1 temp |awk '{a=$(NF-1);b=$NF;if (a>2*b&&a>10) print $0}' >templ
#paste $1 temp |awk '{a=$NF;b=$(NF-1);if (a>2*b&&a>10) print $0}' >tempr
#wc -l templ
#wc -l tempr
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapMinusRawRep1.bigWig $1 1 |awk '{if (NR==1) printf "Minus"; printf "\t"$1} END {printf "\n"}' >templ
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapPlusRawRep1.bigWig $1 1 |awk '{if (NR==1) printf "Plus"; printf "\t"$1 } END {printf "\n"}' >tempr
#/home/ruijie/gbdb/hg19/bbi/wgEncodeRikenCageGm12878CellPapMinusSignalRep1.bigWig
#/home/ruijie/gbdb/hg19/bbi/wgEncodeRikenCageGm12878CellPapPlusSignalRep1.bigWig
mrna="/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig"
prna="/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig"
#mrna="/gbdb/hg19/bbi/wgEncodeRikenCageGm12878CellPapMinusSignalRep1.bigWig"
#prna="/gbdb/hg19/bbi/wgEncodeRikenCageGm12878CellPapPlusSignalRep1.bigWig"

bigWigSummaryBatch $mrna $1 1 |awk '{if (NR==1) printf "LMinus"; printf "\t"$1} END {printf "\n"}' >temp.dat
bigWigSummaryBatch $prna $1 1 |awk '{if (NR==1) printf "LPlus"; printf "\t"$1 } END {printf "\n"}' >>temp.dat
bigWigSummaryBatch $mrna $2 1 |awk '{if (NR==1) printf "MMinus"; printf "\t"$1} END {printf "\n"}' >>temp.dat
bigWigSummaryBatch $prna $2 1 |awk '{if (NR==1) printf "MPlus"; printf "\t"$1 } END {printf "\n"}' >>temp.dat
bigWigSummaryBatch $mrna $3 1 |awk '{if (NR==1) printf "RMinus"; printf "\t"$1} END {printf "\n"}' >>temp.dat
bigWigSummaryBatch $prna $3 1 |awk '{if (NR==1) printf "RPlus"; printf "\t"$1 } END {printf "\n"}' >>temp.dat
#cat templ tempr >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R
