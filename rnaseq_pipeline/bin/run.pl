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
        Usage: $0 GSEXXXXX expr.table group1.lst group2.lst

USAGE
print "$usage";
exit(1);
};

my $gse=shift;
my $expr=shift;
my $group1=shift;
my $group2=shift;

## create table degs
unless(-f "01.degs.finished"){
	`perl $Bin/call_degs.pl $expr $group1 $group2`;
	`touch 01.degs.finished`;
}

## correlation
`Rscript $Bin/correlation.heatmap.r expr.table`;
`convert -density 150 expr.table.pdf -quality 90 expr.table.png`;

`Rscript $Bin/correlation.heatmap.r degs.expr`;
`convert -density 150 degs.expr.pdf -quality 90 degs.expr.png`;

## DEG
`head -1 degs.tab > degs.header.txt`;
`grep -v "^GeneID" degs.tab | awk '\$5 >= 0.58' > up.fc.txt`;
`grep -v "^GeneID" degs.tab | awk '\$5 <= -0.58' > down.fc.txt`;

## just degs
`cat up.fc.txt down.fc.txt > selected.degs.lst`;
`head -1 degs.expr > degs.expr.header`;
`perl $Bin/fishInWinter.pl selected.degs.lst degs.expr > selected.degs.expr.raw`;
`cat degs.expr.header selected.degs.expr.raw > selected.degs.expr`;
`Rscript $Bin/correlation.heatmap.r selected.degs.expr`;
`convert -density 150 selected.degs.expr.pdf -quality 90 selected.degs.expr.png`;

my $mir_count = `cat up.fc.txt down.fc.txt | egrep "(mir|Mir|MIR)" | wc -l`;
chomp($mir_count);

if($mir_count < 30){
	`cut -f 1 up.fc.txt | sort | uniq > up.lst`;
	`cut -f 1 down.fc.txt | sort | uniq > down.lst`;
}elsif($mir_count > 30){
	`perl $Bin/get_MIR_gene.pl up.fc.txt > up.lst`;
	`perl $Bin/get_MIR_gene.pl down.fc.txt > down.lst`;
}
# get the data

my %g1;
my %g2;
my @color;
open IN,"$group1" || die $!;
while(<IN>){
	push @color,"black";
	chomp;
	$g1{$_} = "black";
}
close IN;

open IN,"$group2" || die $!;
while(<IN>){
	push @color,"red";
	chomp;
	$g2{$_} = "red";
}
close IN;

my $color_str = join ",",@color;

open IN,"degs.expr" || die $!;
my $header = <IN>;
chomp($header);
my @header = split /\t/, $header;
my @gi;
my @new_header;
for(my $i=0; $i<scalar(@header); $i++){
	if(exists $g1{$header[$i]}){
		push @gi,$g1{$header[$i]};
		push @new_header,$header[$i];
	}elsif(exists $g2{$header[$i]}){
		push @gi,$g2{$header[$i]};
		push @new_header,$header[$i];
	}
}
close IN;

my $col_str = join ",",@gi;
my $header_str = join ",",@new_header;

print STDERR "$col_str\n";

#print "$col_str\n";
#print "$header_str\n";
### PCA analysis
`Rscript $Bin/pca.R $col_str`;
`convert -density 150 pcaBarplot.pdf -quality 90 pcaBarplot.png`;
`convert -density 150 PCA2d.pdf -quality 90 PCA2d.png`;

print STDERR "$col_str finihsed\n";
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

### degs
#`cut -f 1,2,3,6,7,8 degs.tab > DEGs_table.lst`;
`cp degs.tab DEGs_table.lst`;
`head -1 DEGs_table.lst > DEGs_table.header.txt`;
`perl $Bin/CreatExcelFile.py $Bin/DEGs_table_meta.txt $Bin/header.txt ./DEGs_table.lst ./DEGs_table.xlsx`;

`cat DEGs_table.header.txt up.fc.txt > DEGs_table_up.lst`;
`cat DEGs_table.header.txt down.fc.txt > DEGs_table_down.lst`;

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
			<th>转录本数</th>
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

