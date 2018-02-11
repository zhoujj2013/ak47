# lncRNA
mkdir lncRNA
cd lncRNA
grep "lncRNA" /disk4/database/hg19/refGene.anno | cut -f 3 | sed 's/,/\n/g' | sort | uniq > lncRNA.anno
~/software/cufflinks-2.0.2.Linux_x86_64/gffread /disk4/database/hg19/refGene.gtf -g /disk4/database/hg19/genome.fa -w refgene.fa
perl ~/bin/fishInWinter.pl -bc 1 -ff fasta lncRNA.anno refgene.fa > lncRNA.fa

perl ~/bin/fastaDeal.pl --cutf 8 lncRNA.fa

ls lncRNA.fa.cut/* | while read line;do echo "~/github_project/ak47/S001/bin/miRanda/bin/miranda ../mature.miRNA.fa $line > $line.miranda.output";done > miranda.sh;

perl ~/bin/multi-process.pl -cpu 8 miranda.sh > log 2> err

cat lncRNA.fa.cut/lncRNA.fa.*.miranda.output > lncRNA.fa.miranda.output

egrep "^>\b" ./lncRNA.fa.miranda.output | sed "s/>//g" > ./lncRNA.fa.miranda.result
cd -

