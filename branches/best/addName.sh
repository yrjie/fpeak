cp $1 temp1
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,".",$4}' temp1 >$1
