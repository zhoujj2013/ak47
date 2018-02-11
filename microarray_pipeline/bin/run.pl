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

        Run microarray analysis.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 GSEXXXX GPLXXX 1110000XXXXX

USAGE
print "$usage";
exit(1);
};

my $gse=shift;
my $gpl=shift;
my $sample_string=shift;

unless(-f "01.preprosessing.finished"){
	`Rscript $Bin/preprosessing.R $gse $gpl $sample_string`;
	`touch 01.preprosessing.finished`;	
}

## reformat
`cp deg.tab deg.tab.raw`;
`perl $Bin/replace_trans_id.pl deg.tab.raw $Bin/../data/hg19.refGene.anno > deg.tab`;

## add probe id
open OUT,">","./expr.table" || die $!;
open IN,"./expr.table.raw" || die $!;
my $header = <IN>;
my @header = split /\t/,$header;
unshift @header,"ProbeID";
print OUT join "\t",@header;
while(<IN>){
	print OUT "$_";
}
close IN;
close OUT;

## correlation
`Rscript $Bin/correlation.heatmap.r expr.table`;
`convert -density 150 expr.table.pdf -quality 90 expr.table.png`;

## DEG
`perl $Bin/get_degs.pl deg.tab expr.table 0.58 > degs.expr`;
`cut -f 7 up.fc.txt | sort | uniq | perl -ne 'chomp; my \@t=split \/\\/\/; print join "\\n",\@t;print "\\n";' | sed '/^\$/d' | sort | uniq > up.lst`;
`cut -f 7 down.fc.txt | sort | uniq | perl -ne 'chomp; my \@t=split \/\\/\/; print join "\\n",\@t;print "\\n";' | sed '/^\$/d' | sort | uniq > down.lst`;

`Rscript $Bin/correlation.heatmap.r degs.expr`;
`convert -density 150 degs.expr.pdf -quality 90 degs.expr.png`;

# degs
`cat up.fc.txt down.fc.txt > selected.degs.lst`;
`head -1 degs.expr > degs.expr.header`;
`perl $Bin/fishInWinter.pl selected.degs.lst degs.expr > selected.degs.expr.raw`;
`cat degs.expr.header selected.degs.expr.raw > selected.degs.expr`;
`Rscript $Bin/correlation.heatmap.r selected.degs.expr`;
`convert -density 150 selected.degs.expr.pdf -quality 90 selected.degs.expr.png`;

# get the data
my @markList = split //,$sample_string;
my %s;
my $i=0;
open IN,"./sample.lst" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$s{$t[0]} = $markList[$i];
	$i++;
}
close IN;

open IN,"degs.expr" || die $!;
$header = <IN>;
chomp($header);
@header = split /\t/,$header;
shift @header;
close IN;
my $header_str = join ",",@header;


my @col;
my $j = 1;
foreach my $c (@header){
	$j++;
	push @col,$j;
}
my $col_str = join ",",@col;

my @color;
foreach my $h (@header){
	if($s{$h} eq "0"){
		push @color,"black";
	}elsif($s{$h} eq "1"){
		push @color,"red";
	}
}

my $color_str = join ",",@color;

## PCA analysis
`Rscript $Bin/pca.R $color_str`;
`convert -density 150 pcaBarplot.pdf -quality 90 pcaBarplot.png`;
`convert -density 150 PCA2d.pdf -quality 90 PCA2d.png`;


# annotation
unless(-f "02.anno.finished"){
	`sh ~/github_project/jjUtil/DavidAnno/anno.sh up.lst hsapiens up`;
	`sh ~/github_project/jjUtil/DavidAnno/anno.sh down.lst hsapiens down`;
	`touch 02.anno.finished`;
}

## deg clustering
#print "Rscript $Bin/drawHeatmap4AllSamplesColRow.r cluster $col_str $header_str degs.expr\n";
`Rscript $Bin/drawHeatmap4AllSamplesColRow.r cluster $col_str $header_str degs.expr`;
`convert -density 150 cluster.all.sample.heatmap.pdf -quality 90 cluster.all.sample.heatmap.png`;

## create excel file for report
### expression
`perl $Bin/CreatExcelFile.py $Bin/samples_anno_meta.txt $Bin/header.txt ./samples_anno.lst ./samples_anno.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/Gene_expression_table_meta.txt $Bin/header.txt ./expr.table ./Gene_expression_table.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/Gene_expression_table_meta.txt $Bin/header.txt ./degs.expr ./DEGs_expr_table.xlsx`;

## degs
`cut -f 1,2,3,6,7,8 deg.tab > DEGs_table.lst`;
`head -1 DEGs_table.lst > DEGs_table.header.txt`;
`perl $Bin/CreatExcelFile.py $Bin/DEGs_table_meta.txt $Bin/header.txt ./DEGs_table.lst ./DEGs_table.xlsx`;

`cut -f 1,2,3,6,7,8 up.fc.txt | cat DEGs_table.header.txt - > DEGs_table_up.lst`;
`cut -f 1,2,3,6,7,8 down.fc.txt | cat DEGs_table.header.txt - > DEGs_table_down.lst`;

`perl $Bin/CreatExcelFile.py $Bin/DEGs_table_meta.txt $Bin/header.txt ./DEGs_table_up.lst ./DEGs_table_up.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/DEGs_table_meta.txt $Bin/header.txt ./DEGs_table_down.lst ./DEGs_table_down.xlsx`;

### annotation
`head -1 down.david/down.david/down.input.tableReport.final.txt > anno.header.txt`;
`less -S down.david/down.david/down.input.tableReport.final.txt | grep KEGG | cat anno.header.txt - > DEGs_anno_kegg_down.lst`;
`less -S up.david/up.david/up.input.tableReport.final.txt | grep KEGG | cat anno.header.txt - > DEGs_anno_kegg_up.lst`;

`less -S down.david/down.david/down.input.tableReport.final.txt | grep GOTERM | cat anno.header.txt - > DEGs_anno_GO_down.lst`;
`less -S up.david/up.david/up.input.tableReport.final.txt | grep GOTERM | cat anno.header.txt - > DEGs_anno_GO_up.lst`;

`perl $Bin/CreatExcelFile.py $Bin/anno.meta.txt $Bin/header.txt DEGs_anno_kegg_down.lst ./DEGs_anno_kegg_down.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/anno.meta.txt $Bin/header.txt DEGs_anno_kegg_up.lst ./DEGs_anno_kegg_up.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/anno.meta.txt $Bin/header.txt DEGs_anno_GO_down.lst ./DEGs_anno_GO_down.xlsx`;
`perl $Bin/CreatExcelFile.py $Bin/anno.meta.txt $Bin/header.txt DEGs_anno_GO_up.lst ./DEGs_anno_GO_up.xlsx`;

# create GO and KEGG graph
`grep GOTERM up.david/up.david/up.input.tableReport.final.txt | cut -f 4,12 | awk -F"\\t" '{print \$2"\\t"\$1}' | head -10  > up.GO.txt`;
`Rscript $Bin/anno.barplot.r up.GO.txt up.GO.pdf`;
`convert -density 150 up.GO.pdf -quality 90 up.GO.png`;

`grep GOTERM down.david/down.david/down.input.tableReport.final.txt | cut -f 4,12 | awk -F"\\t" '{print \$2"\\t"\$1}' | head -10  > down.GO.txt`;
`Rscript $Bin/anno.barplot.r down.GO.txt down.GO.pdf`;
`convert -density 150 down.GO.pdf -quality 90 down.GO.png`;

`grep KEGG_PATHWAY up.david/up.david/up.input.tableReport.final.txt | cut -f 4,12 | awk -F"\\t" '{print \$2"\\t"\$1}' | head -10  > up.KEGG.txt`;
`Rscript $Bin/anno.barplot.r up.KEGG.txt up.KEGG.pdf`;
`convert -density 150 up.KEGG.pdf -quality 90 up.KEGG.png`;

`grep KEGG_PATHWAY down.david/down.david/down.input.tableReport.final.txt | cut -f 4,12 | awk -F"\\t" '{print \$2"\\t"\$1}' | head -10  > down.KEGG.txt`;
`Rscript $Bin/anno.barplot.r down.KEGG.txt down.KEGG.pdf`;
`convert -density 150 down.KEGG.pdf -quality 90 down.KEGG.png`;

# create link file
`perl $Bin/kegg_link/getPATHWAYmap.pl up.david/up.david/up.input.tableReport.final.txt down.david/down.david/down.input.tableReport.final.txt $Bin/kegg_link/human.geneid.lst > kegg.link.txt`;

### to arrange the result
`rm -r report` if(-d "report");
mkdir "report" unless(-d "report");
mkdir "report/result" unless(-d "report/result");
open IN,"$Bin/result.txt" || die $!;
while(<IN>){
	chomp;
	`cp $_ report/result/`;
}
close IN;

`cp -r $Bin/report_html/* report/`;

### to create html file
`cp report/02.html report/02.html.bk`;
open IN,"report/02.html.bk" || die $!;
open OUT,">","report/02.html" || die $!;
while(<IN>){
	my $line = $_;
	if(/%GSEID%/){
		$line =~ s/%GSEID%/$gse/g;
		print OUT "$line";
	}else{
		print OUT "$line";
	}
}
close IN;
close OUT;

# 03.html
my $up_probe_num = `wc -l up.fc.txt | awk '{print \$1}'`;
chomp($up_probe_num);
my $down_probe_num = `wc -l down.fc.txt | awk '{print \$1}'`;
chomp($down_probe_num);

my $up_gene_num = `wc -l up.lst | awk '{print \$1}'`;
chomp($up_gene_num);
my $down_gene_num = `wc -l down.lst | awk '{print \$1}'`;
chomp($down_gene_num);

my $table = << "TABLE";
<div>
<table class="table table-hover">
	<thead>
		<tr>
			<th>上调/下调</th>
			<th>探针数</th>
			<th>基因数</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td>上调（up）</td>
			<td>$up_probe_num</td>
			<td>$up_gene_num</td>
		</tr>
		<tr>
			<td>下调（down）</td>
			<td>$down_probe_num</td>
			<td>$down_gene_num</td>
		</tr>
	</tbody>
</table>
</div>
<p></p>
TABLE

`cp report/03.html report/03.html.bk`;
open IN,"report/03.html.bk" || die $!;
open OUT,">","report/03.html" || die $!;
while(<IN>){
    my $line = $_;
    if(/%EXPRTABLE%/){
        $line =~ s/%EXPRTABLE%/$table/g;
        print OUT "$line";
    }else{
        print OUT "$line";
    }
}
close IN;
close OUT;

# 05.html

my %klink;
open IN,"./kegg.link.txt" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$klink{$t[0]} = $t[1];
}
close IN;

my $uplink ="";
open IN,"./up.KEGG.txt" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $outline = "<p><a href=\"$klink{$t[0]}\">$t[0]</a></p>";
	$uplink .= "$outline\n";
}
close IN;

my $downlink ="";
open IN,"./down.KEGG.txt" || die $!;
while(<IN>){
	chomp;
    my @t = split /\t/;
    my $outline = "<p><a href=\"$klink{$t[0]}\">$t[0]</a></p>";
    $downlink .= "$outline\n";
	#print STDERR "$outline\n";
}
close IN;

`cp report/05.html report/05.html.bk`;
open IN,"report/05.html.bk" || die $!;
open OUT,">","report/05.html" || die $!;
while(<IN>){
    my $line = $_;
    if(/%UPLINK%/){
        $line =~ s/%UPLINK%/$uplink/g;
        print OUT "$line";
    }elsif(/%DOWNLINK%/){
		$line =~ s/%DOWNLINK%/$downlink/g;
		print OUT "$line";
	}else{
        print OUT "$line";
    }
}
close IN;
close OUT;

`rm report/03.html.bk report/02.html.bk report/05.html.bk`;

