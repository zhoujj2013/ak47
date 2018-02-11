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

my $f=shift;

my %utr;
open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	push @{$utr{$t[3]}{utr}},\@t;
	$utr{$t[3]}{strand} = $t[5];
}
close IN;

foreach my $trans_id (keys %utr){
	if($utr{$trans_id}{strand} eq "-"){
		my @utr_sorted = sort {$a->[1] <=> $b->[1]} @{$utr{$trans_id}{utr}};
		my @new_utr;
		my $j = scalar(@utr_sorted)-1;
		for(my $i = scalar(@utr_sorted)-1; $i >= 0; $i--){
			my $h = $j - $i;
			$utr_sorted[$i][3] = "$utr_sorted[$i][3]"."-$h";
			#push @new_utr,$utr_sorted[$i];
			print join "\t",@{$utr_sorted[$i]};
			print "\n";
		}
	}elsif($utr{$trans_id}{strand} eq "+"){	
		my @utr_sorted = sort {$a->[1] <=> $b->[1]} @{$utr{$trans_id}{utr}};
		my @new_utr;
		#my $j = scalar(@utr_sorted)-1;
		for(my $i = 0; $i < scalar(@utr_sorted); $i++){
			#my $h = $j - $i;
			$utr_sorted[$i][3] = "$utr_sorted[$i][3]"."-$i";
			#push @new_utr,$utr_sorted[$i];
			print join "\t",@{$utr_sorted[$i]};
			print "\n";
		}
	}
}

