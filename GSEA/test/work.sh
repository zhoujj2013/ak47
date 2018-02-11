/disk4/project/P20170819/GSE10616/ColonOnly_CD.vs.healthy/deg.tab.raw | cut -f 6,7 | awk '{print $2"\t"$1}' | grep -v "Gene" | sort -k2gr > rank.lst

sh ../bin/run.GSEA.human.sh <expr.foldchange.rnk> <prefix>


