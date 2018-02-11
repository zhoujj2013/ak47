# utr3.bed from buildDB pipeline
mkdir mRNA
cd mRNA
perl $(dirname $0)/Rename3UTR.pl ../utr3.bed > utr3.renamed.bed
bedtools getfasta -fi /disk4/database/hg19/genome.fa -bed ./utr3.renamed.bed -fo ./utr3.renamed.fa -s -name 2> err
perl $(dirname $0)/link3UTR.pl utr3.renamed.fa > utr3.renamed.link.fa

ln -s utr3.renamed.link.fa utr3.fa

perl ~/bin/fastaDeal.pl --cutf 8 utr3.fa

ls utr3.fa.cut/* | while read line;do echo "~/github_project/ak47/S001/bin/miRanda/bin/miranda ../mature.miRNA.fa $line > $line.miranda.output";done > miranda.sh;

perl ~/bin/multi-process.pl -cpu 8 miranda.sh > log 2> err

cat utr3.fa.cut/utr3.fa.*.miranda.output > utr3.fa.miranda.output

rm -r utr3.fa.cut

egrep "^>\b" ./utr3.fa.miranda.output | sed "s/>//g" > ./utr3.fa.miranda.result
cd -

