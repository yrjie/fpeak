# $1 plusSig
# $2 minusSig
# $3 prefix
if [ $# -lt 3 ]
then
    echo "Usage: sh callPeak.sh plusSig minusSig prefix"
    exit
fi
plusSig=$1
minusSig=$2
prefix=$3
plusBed=$prefix'Plus.bed'
minusBed=$prefix'Minus.bed'
peakBed=$prefix'Peak.bed'
top=200000
sort -k5 -gr $plusSig |head -n $top >$plusBed
sort -k5 -gr $minusSig |head -n $top >$minusBed
cat $plusBed $minusBed |sortBed |mergeBed -scores mean >$peakBed
