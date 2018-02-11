#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
use List::Util qw(sum);
use Statistics::R;

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script create makefile for LncFunNet analysis.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 expr.table g1 g2
        
USAGE
print "$usage";
exit(1);
};

my $expr_f=shift;
my $group1 = shift;
my $group2 = shift;

my %g1;
open IN,"$group1" || die $!;
while(<IN>){
    chomp;
    my @t = split /\t/;
    $g1{$t[0]} = 1;
}
close IN;

my %g2;
open IN,"$group2" || die $!;
while(<IN>){
    chomp;
    my @t = split /\t/;
    $g2{$t[0]} = 1;
}
close IN;

open IN,"$expr_f" || die $!;
open OUT,">","degs.expr" || die $!;
open OUT1,">","degs.tab" || die $!;
my $header = <IN>;
chomp($header);
my @header = split /\t/,$header;
my $first = $header[0];
my %gg1;
my %gg2;
my @reheader;
for(my $i=0; $i < scalar(@header); $i++){
	if(exists $g1{$header[$i]}){
		$gg1{$i} = $header[$i];
		push @reheader,$header[$i];
	}
	if(exists $g2{$header[$i]}){
        $gg2{$i} = $header[$i];
		push @reheader,$header[$i];
    }
}

print OUT "$first\t";
print OUT join "\t",@reheader;
print OUT "\n";

print OUT1 "GeneID\tgroup1_mean\tgroup2_mean\tFC\tlogFC\tpvalue\tadj_pvalue\n";
my @out1;
my @pvalue;

open G1OUT,">","g1.expr" || die $!;
open G2OUT,">","g2.expr" || die $!;

my @geneid;

while(<IN>){
	chomp;
	my @t = split /\t/;
	push @geneid,$t[0];
	my @g1;
	my @g2;
	for(my $j = 0; $j < scalar(@t); $j++){
		if(exists $gg1{$j}){
			push @g1,$t[$j];
		}
		if(exists $gg2{$j}){
            push @g2,$t[$j];
        }
	}
	print G1OUT join "\t",@g1;
	print G1OUT "\n";
	print G2OUT join "\t",@g2;
	print G2OUT "\n";
	
	print OUT "$t[0]\t";
	print OUT join "\t",@g1;
	print OUT "\t";
	print OUT join "\t",@g2;
	print OUT "\n";
}
close IN;
close OUT;

`Rscript $Bin/ttest.r`;

my $i = 0;
open IN,"./ttest.result" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	unshift @t,$geneid[$i];
	print OUT1 join "\t",@t;
	print OUT1 "\n";
	$i++;
}
close IN;
close OUT1;


#my $R = Statistics::R->new();
#my $str = join ",",@pvalue;
#$R->run(qq`adjp=p.adjust(c($str), method="hochberg")`);
#my $adj_pvalue = $R->get('adjp');
#print Dumper($adj_pvalue);
#$R->stop();
#undef $R;

#my $k = 0;
##print STDERR $adj_pvalue->[0],"\n";
#foreach my $o (@out1){
#	#print STDERR "$adj_pvalue->[$k]\n";
#	push @{$o},$adj_pvalue->[$k];
#	print OUT1 join "\t",@{$o};
#	print OUT1 "\n";
#	$k++;
#}
#
#close OUT1;
#
#sub ttest{
#	my $ggg1=shift;
#	my $ggg2=shift;
#	
#	my $R = Statistics::R->new();
#	
#	my $a=join ",",@{$ggg1};
#	my $b=join ",",@{$ggg2};
#	$R->run(qq`r=t.test(c($a), c($b),alternative = "two.sided")`);
#	$R->run(q`y <- r$p.value`);
#	my $p_value = $R->get('y');
#	$R->stop();
#	undef $R;
#	return $p_value
#}
