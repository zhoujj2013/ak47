perl ../bin/run_alignment.pl config.txt sample.lst
OUTDIR  ./out
SAMPLE  ./sample.lst
thread  4
aligner tophat2
tophat2 /home/zhoujj/software/tophat-2.1.1.Linux_x86_64/tophat2
bowtie2 /usr/bin/bowtie2
bwa     /home/zhoujj/software/bwa-0.6.2/bwa
soap    /home/zhoujj/software/soap2.21release/soap
samtools        /usr/bin/samtools
picard  /home/zhoujj/software/picard-tools-2.2.1/picard.jar
star    /home/zhoujj/software/STAR-2.5.2b/bin/Linux_x86_64_static/STAR
gtf     /disk4/project/LycPUFA/03CheckTranscriptome/spe/mm10.gencode.vM14.annotation.gtf
index   /disk4/project/LycPUFA/03CheckTranscriptome/spe/mm10
