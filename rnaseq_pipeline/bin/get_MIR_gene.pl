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

my %gene;
open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my @mir = split /-/,$t[0];
	my $id = "";
	if($t[0] =~ /let/){
		$id="MIRLET$mir[2]";
	}else{
		$id="MIR$mir[2]";
	}
	$gene{$id}=1;
}
close IN;


foreach my $k (keys %gene){
	$k=uc($k);
	print "$k\n";
}
