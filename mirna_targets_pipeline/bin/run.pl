#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        Run microRNA target pipeline.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 up_mir_f down_mir_f up_lncrna_f down_lncrna_f up_mrna_f down_mrna_f mir_mrna_int_f mir_lncrna_int_f mir_mrna_BS_f mir_lncrna_BS_f mirTarbase_f

USAGE
print "$usage";
exit(1);
};

my $up_mir_f=shift;
my $down_mir_f=shift;

my $up_lncrna_f = shift;
my $down_lncrna_f = shift;

my $up_mrna_f = shift;
my $down_mrna_f = shift;

my $mir_mrna_f = shift;
my $mir_lncrna_f = shift;

my $mir_mrna_bs_f = shift;
my $mir_lncrna_bs_f = shift;

my $mirTarbase_f = shift;

# for lncRNA
`perl ~/bin/fishInWinter.pl $up_mir_f  $mir_lncrna_f | perl ~/bin/fishInWinter.pl -fc 2 $down_lncrna_f - > up.mirna.vs.down.lncrna.cl.int`;
`cut -f 1,2 up.mirna.vs.down.lncrna.cl.int | sort | uniq -c | awk '{print \$2"\t"\$3"\t"\$1}' > up.mirna.vs.down.lncrna.cl.int.rmdup`;

`perl ~/bin/fishInWinter.pl $down_mir_f $mir_lncrna_f | perl ~/bin/fishInWinter.pl -fc 2 $up_lncrna_f - > down.mirna.vs.up.lncrna.cl.int`;
`cut -f 1,2 down.mirna.vs.up.lncrna.cl.int | sort | uniq -c | awk '{print \$2"\t"\$3"\t"\$1}' > down.mirna.vs.up.lncrna.cl.int.rmdup`;

# for mRNA
`perl ~/bin/fishInWinter.pl $up_mir_f  $mir_mrna_f | perl ~/bin/fishInWinter.pl -fc 2 $down_mrna_f - > up.mirna.vs.down.mrna.cl.int`;
`cut -f 1,2 up.mirna.vs.down.mrna.cl.int | sort | uniq -c | awk '{print \$2"\t"\$3"\t"\$1}' > up.mirna.vs.down.mrna.cl.int.rmdup`;

`perl ~/bin/fishInWinter.pl $down_mir_f $mir_mrna_f | perl ~/bin/fishInWinter.pl -fc 2 $up_mrna_f - > down.mirna.vs.up.mrna.cl.int`;
`cut -f 1,2 down.mirna.vs.up.mrna.cl.int | sort | uniq -c | awk '{print \$2"\t"\$3"\t"\$1}' > down.mirna.vs.up.mrna.cl.int.rmdup`;

# combine information mRNA
`perl $Bin/combine_information_mrna.pl up.mirna.vs.down.mrna.cl.int.rmdup $down_mrna_f $mirTarbase_f > up.mirna.vs.down.mrna.cl.int.rmdup.txt 2>up.mirna.vs.down.mrna.cl.int.rmdup.err`;

`perl $Bin/combine_information_mrna.pl down.mirna.vs.up.mrna.cl.int.rmdup $up_mrna_f $mirTarbase_f > down.mirna.vs.up.mrna.cl.int.rmdup.txt 2>down.mirna.vs.up.mrna.cl.int.rmdup.err`;

# combine information lncRNA
`perl $Bin/combine_information_lncrna.pl up.mirna.vs.down.lncrna.cl.int.rmdup $down_lncrna_f $mirTarbase_f > up.mirna.vs.down.lncrna.cl.int.rmdup.txt 2>up.mirna.vs.down.lncrna.cl.int.rmdup.err`;

`perl $Bin/combine_information_lncrna.pl down.mirna.vs.up.lncrna.cl.int.rmdup $up_lncrna_f $mirTarbase_f > down.mirna.vs.up.lncrna.cl.int.rmdup.txt 2>down.mirna.vs.up.lncrna.cl.int.rmdup.err`;

# for input file to excel
`cat $Bin/mir.input.header $up_mir_f > Up_regulated.mir.lst`;
`cat $Bin/mir.input.header $down_mir_f > Down_regulated.mir.lst`;

`cat $Bin/gene.input.header $up_lncrna_f > Up_regulated.lncrna.lst`;
`cat $Bin/gene.input.header $down_lncrna_f > Down_regulated.lncrna.lst`;

`cat $Bin/gene.input.header $up_mrna_f > Up_regulated.mrna.lst`;
`cat $Bin/gene.input.header $down_mrna_f > Down_regulated.mrna.lst`;

`perl $Bin/CreatExcelFile.py $Bin/mir.DEGs_table_meta.txt $Bin/header.txt Up_regulated.mir.lst Up_regulated.mir.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/mir.DEGs_table_meta.txt $Bin/header.txt Down_regulated.mir.lst Down_regulated.mir.xlsx`;

`perl $Bin/CreatExcelFile.py $Bin/gene.DEGs_table_meta.txt $Bin/header.txt Up_regulated.lncrna.lst Up_regulated.lncrna.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/gene.DEGs_table_meta.txt $Bin/header.txt Down_regulated.lncrna.lst Down_regulated.lncrna.xlsx`;

`perl $Bin/CreatExcelFile.py $Bin/gene.DEGs_table_meta.txt $Bin/header.txt Up_regulated.mrna.lst Up_regulated.mrna.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/gene.DEGs_table_meta.txt $Bin/header.txt Down_regulated.mrna.lst Down_regulated.mrna.xlsx`;

# for result

`cat $Bin/result.header down.mirna.vs.up.lncrna.cl.int.rmdup.txt > down.mirna.vs.up.lncrna.txt`;
`cat $Bin/result.header down.mirna.vs.up.mrna.cl.int.rmdup.txt > down.mirna.vs.up.mrna.txt`;
`cat $Bin/result.header up.mirna.vs.down.lncrna.cl.int.rmdup.txt > up.mirna.vs.down.lncrna.txt`;
`cat $Bin/result.header up.mirna.vs.down.mrna.cl.int.rmdup.txt > up.mirna.vs.down.mrna.txt`;

`perl $Bin/CreatExcelFile.py $Bin/result.meta $Bin/header.txt down.mirna.vs.up.lncrna.txt down.mirna.vs.up.lncrna.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/result.meta $Bin/header.txt down.mirna.vs.up.mrna.txt down.mirna.vs.up.mrna.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/result.meta $Bin/header.txt up.mirna.vs.down.lncrna.txt  up.mirna.vs.down.lncrna.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/result.meta $Bin/header.txt up.mirna.vs.down.mrna.txt  up.mirna.vs.down.mrna.xlsx`;

# get binding site
`perl $Bin/get_binding_site.pl down.mirna.vs.up.mrna.cl.int.rmdup.txt $mir_mrna_bs_f > down.mirna.vs.up.mrna.bindingSite.txt`;
`perl $Bin/get_binding_site.pl down.mirna.vs.up.lncrna.cl.int.rmdup.txt $mir_lncrna_bs_f > down.mirna.vs.up.lncrna.bindingSite.txt`;
`perl $Bin/get_binding_site.pl up.mirna.vs.down.lncrna.cl.int.rmdup.txt $mir_lncrna_bs_f > up.mirna.vs.down.lncrna.bindingSite.txt`;
`perl $Bin/get_binding_site.pl up.mirna.vs.down.mrna.cl.int.rmdup.txt $mir_mrna_bs_f > up.mirna.vs.down.mrna.bindingSite.txt`;

# create a network
`perl $Bin/create_igraph_rscript.pl down.mirna.vs.up.mrna.cl.int.rmdup.txt 10 down.mirna.vs.up.mrna > down.mirna.vs.up.mrna.R`;
`Rscript down.mirna.vs.up.mrna.R`;
`convert -density 150 down.mirna.vs.up.mrna.pdf -quality 90 down.mirna.vs.up.mrna.png`;

`perl $Bin/create_igraph_rscript.pl down.mirna.vs.up.lncrna.cl.int.rmdup.txt 10 down.mirna.vs.up.lncrna > down.mirna.vs.up.lncrna.R`;
`Rscript down.mirna.vs.up.lncrna.R`;
`convert -density 150 down.mirna.vs.up.lncrna.pdf  -quality 90 down.mirna.vs.up.lncrna.png`;

`perl $Bin/create_igraph_rscript.pl up.mirna.vs.down.lncrna.cl.int.rmdup.txt 10 up.mirna.vs.down.lncrna > up.mirna.vs.down.lncrna.R`;
`Rscript up.mirna.vs.down.lncrna.R`;
`convert -density 150 up.mirna.vs.down.lncrna.pdf -quality 90 up.mirna.vs.down.lncrna.png`;

`perl $Bin/create_igraph_rscript.pl up.mirna.vs.down.mrna.cl.int.rmdup.txt  10 up.mirna.vs.down.mrna > up.mirna.vs.down.mrna.R`;
`Rscript up.mirna.vs.down.mrna.R`;
`convert -density 150 up.mirna.vs.down.mrna.pdf -quality 90 up.mirna.vs.down.mrna.png`;

# create report
`rm -r ./report` if(-d "./report");
`cp -r $Bin/report_html ./report`;
mkdir "./report/result" unless(-d "./report/result");
open IN,"$Bin/result.txt" || die $!;
while(<IN>){
	chomp;
	`cp $_ ./report/result/`;
}
close IN;

