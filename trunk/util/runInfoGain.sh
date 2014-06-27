awk 'BEGIN{OFS="\t"}{if (NR==1) next; for (i=3;i<=7;i++) $i/=$8; print $3,$4,$5,$6,$7,$8}' uncovTFMerg_frag.xls >pos.tab
awk 'BEGIN{OFS="\t"}{if (NR==1) next; for (i=3;i<=7;i++) $i/=$8; print $3,$4,$5,$6,$7,$8}' AtacGm12878_nonShortNoTFTop_frag.xls >neg.tab

python ~/bin/tab2arff.py pos.tab neg.tab nonShort >temp.arff
sh ~/bin/runSelector.sh

sh ~/bin/runSelector.sh -v temp.arff
