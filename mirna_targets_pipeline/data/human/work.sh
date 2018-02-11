ln -s $lncfuntk_db/miRNAdb/mature.miRNA.fa
ln -s $lncfuntk_db/utr3.bed

sh $PWD/bin/mrna_run.sh
sh $PWD/bin/lncrna_run.sh
