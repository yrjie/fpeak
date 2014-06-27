#cat AtacGm12878_*|windowBed -a AtacGm12878Rmdup100NoStL100Ex.bed -b stdin -v -w 0 >unique100.bed
#cat AtacGm12878_150-200NoStL100Ex.bed AtacGm12878_200-300NoStL100Ex.bed AtacGm12878_300+NoStL100Ex.bed AtacGm12878Rmdup100NoStL100Ex.bed |windowBed -a AtacGm12878_100-150NoStL100Ex.bed -b stdin -v -w 0 >unique100-150.bed
cat AtacGm12878_100-150NoStL100Ex.bed AtacGm12878_200-300NoStL100Ex.bed AtacGm12878_300+NoStL100Ex.bed AtacGm12878Rmdup100NoStL100Ex.bed |windowBed -a AtacGm12878_150-200NoStL100Ex.bed -b stdin -v -w 0 >unique150-200.bed
cat AtacGm12878_100-150NoStL100Ex.bed AtacGm12878_150-200NoStL100Ex.bed AtacGm12878_300+NoStL100Ex.bed AtacGm12878Rmdup100NoStL100Ex.bed |windowBed -a AtacGm12878_200-300NoStL100Ex.bed -b stdin -v -w 0 >unique200-300.bed
cat AtacGm12878_100-150NoStL100Ex.bed AtacGm12878_150-200NoStL100Ex.bed AtacGm12878_200-300NoStL100Ex.bed AtacGm12878Rmdup100NoStL100Ex.bed |windowBed -a AtacGm12878_300+NoStL100Ex.bed -b stdin -v -w 0 >unique300+.bed
