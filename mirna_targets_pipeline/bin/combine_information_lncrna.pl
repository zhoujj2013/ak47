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
        Usage: $0 int_f lncrna_f tarbase
        
        NOTE: please check config.cfg format in ./demo directory.

USAGE
print "$usage";
exit(1);
};

my $int_f = shift;
my $lncrna_f = shift;
my $tarbase = shift;

my %coding;
open IN,"$lncrna_f" || die $!;

my $header1 = <IN>;
chomp($header1);

while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[0];
	$coding{$id} = \@t;
}
close IN;

my %int;
open IN,"$tarbase" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$int{$t[1]}{$t[3]} = $t[-1];
}
close IN;

#print Dumper(\%int);
open IN,"$int_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $tran_id = $t[1];
	if(exists $coding{$tran_id}){
		my $l = $coding{$tran_id};
		my $sb = $l->[0];
		if(exists $int{$t[0]}{$sb}){
			my $link = "https://www.ncbi.nlm.nih.gov/pubmed/$int{$t[0]}{$sb}";
			print join "\t",@t;
			print "\t";
			my @tmp = @{$coding{$tran_id}};
			shift @tmp;
			print join "\t",@tmp;
			print "\t";
			print "$int{$t[0]}{$sb}\t$link\n";
		}else{	
			print join "\t",@t;
			print "\t";
			my @tmp = @{$coding{$tran_id}};
			shift @tmp;
			print join "\t",@tmp;
			print "\t";
			print "NA\tNA\n";
		}
	}else{
		print STDERR "Missing $tran_id\n";
	}
}
close IN;

