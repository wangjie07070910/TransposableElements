paste <(awk -vOFS="\t" '{print $1,$2,$3,$4}' ${genomecovfilenohead}) *TEsrefined > TEsrefinedNA.genomecovnohead

cat <(ls *TEsrefined | sed 's/\.TEsrefined//g' | awk 'BEGIN{printf "%s\t%s\t%s\t%s\t", "chr","start","end","TE"}{printf "%s\t", $0}END{print "";}') TEsrefinedNA.genomecovnohead  > TEsrefinedNA.genomecov
