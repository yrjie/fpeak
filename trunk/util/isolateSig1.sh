if [ $# -lt 1 ]
then
    echo "Usage: windowSize"
    exit
fi

tmpF=tempIso.bed
echo $1
#awk -v thre=$1 'BEGIN{pre=0;chr="";st=0;end=0;sc=0;OFS="\t"}{if (sc==1&&st-pre>thre&&$2-end>thre) print chr,st,end,sc; pre=end;chr=$1;st=$2;end=$3;sc=$4;}' pilePeak150-_0.bed >$tmpF
awk -v thre=$1 'BEGIN{pre=0;chr="";st=0;end=0;sc=0;OFS="\t"}{if (sc==1&&(st-pre<thre||$2-end<thre)) print chr,st,end,sc; pre=end;chr=$1;st=$2;end=$3;sc=$4;}' pilePeak150-_0.bed >$tmpF
a=`wc -l $tmpF`
b=`windowBed -a $tmpF -b ../allGm12878TFrmdup.bed -u -w 0|wc -l`
echo $a $b

#awk -v thre=$1 'BEGIN{pre=0;chr="";st=0;end=0;sc=0;OFS="\t"}{if (sc==1&&st-pre>thre&&$2-end>thre) print chr,st,end,sc; pre=end;chr=$1;st=$2;end=$3;sc=$4;}' pilePeak150-_0.bed |windowBed -b stdin -a ../allGm12878TFrmdup.bed -u -w 0|wc -l

#awk -v thre=200 'BEGIN{pre=0;chr="";st=0;end=0;sc=0;OFS="\t"}{if (sc==1&&(st-pre<thre||$2-end<thre)) print chr,st,end,sc; pre=end;chr=$1;st=$2;end=$3;sc=$4;}' pilePeak150-_0.bed |wc -l
