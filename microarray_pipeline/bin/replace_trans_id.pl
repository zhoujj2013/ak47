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

        This script create makefile for LncFunNet analysis.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 config.cfg
        
        NOTE: please check config.cfg format in ./demo directory.

USAGE
print "$usage";
exit(1);
};

my $deg_f=shift;
my $refgene_f = shift;

my $trans_count = 0;
my $total = 0;
open IN,"$deg_f" || die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[6];
	next unless(defined $id);
	$trans_count++ if($id =~ /NM/);	
	$total++;
}
close IN;

my %trans2gene;
open IN,"$refgene_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my @tids = split /,/,$t[2];
	foreach my $ttid (@tids){
		$trans2gene{$ttid} = $t[0];
	}
}
close IN;

if($trans_count/$total > 0.5){
	open IN,"$deg_f" || die $!;
	my $header = <IN>;
	chomp($header);
	my @header = split /\t/,$header;
	$header[6] = "Gene.symbol" if($header[6] ne "Gene.symbol");
	$header[7] = "Gene.title" if($header[7] ne "Gene.title");
	print join "\t",@header;
	print "\n";
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $tttid = $t[6];
		if(! exists $trans2gene{$tttid}){
			#print STDERR "$tttid\n";
			$t[7] = "$tttid, $t[7]";
		}else{
			my $gene_symbol = $trans2gene{$tttid};
			$t[6] = $gene_symbol;
			$t[7] = "$tttid, $t[7]";
		}
		print join "\t",@t;
		print "\n";
	}
	close IN;
}else{
	open IN,"$deg_f" || die $!;	
	while(<IN>){
		print "$_";
	}
	close IN;
}
