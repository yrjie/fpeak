awk '{bp=6;for (i=1;i<=bp;i++) printf ($i-($(i+1)+$(i+2))/2)"\t"; for (i=21-bp;i<=20;i++) printf ($i-($(i-1)+$(i-2))/2)"\t"; printf "\n"}' data/Gm12878FaireSp1_c2_ave.sig >templ
awk '{bp=6;for (i=1;i<=bp;i++) printf ($i-($(i+1)+$(i+2))/2)"\t"; for (i=21-bp;i<=20;i++) printf ($i-($(i-1)+$(i-2))/2)"\t"; printf "\n"}' data/Gm12878FaireSp1_c1_ave.sig >tempr
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t";printf "\n"}' data/Gm12878FaireSp1_c2.sig >templ
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t"; printf "\n"}' data/Gm12878FaireSp1_c1.sig >tempr
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t";printf "\n"}' data/Gm12878Sp1By79_c3.sig >templ
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t"; printf "\n"}' data/Gm12878Sp1By79_c4.sig >tempr
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t";printf "\n"}' templ.sig >templ
#awk '{bp=9;for (i=1;i<=bp;i++) printf $i"\t";for (i=21-bp;i<=20;i++) printf $i"\t"; printf "\n"}' tempr.sig >tempr
#awk '{bp=8;for (i=1;i<=bp;i++) printf ($i-($(i+1)+$(i+2))/2)"\t"; for (i=21-bp;i<=20;i++) printf ($i-($(i-1)+$(i-2))/2)"\t"; printf "\n"}' templ.sig >templ
#awk '{bp=8;for (i=1;i<=bp;i++) printf ($i-($(i+1)+$(i+2))/2)"\t"; for (i=21-bp;i<=20;i++) printf ($i-($(i-1)+$(i-2))/2)"\t"; printf "\n"}' tempr.sig >tempr
python ~/bin/plotCorr.py templ
python ~/bin/plotCorr.py tempr
cat templ tempr >temp.sig
python ~/bin/plotCorr.py temp.sig
