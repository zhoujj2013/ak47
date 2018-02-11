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

my %seq;
open IN,"$f" || die $!;
$/ = ">"; <IN>; $/ = "\n";
while(<IN>){
	chomp;
	my $id = $1 if(/(\S+)/);
	my @id = split /-/,$id;
	$/ = ">";
	my $sequence = <IN>;
	chomp($sequence);
	$/ = "\n";
	chomp($sequence);
	$seq{$id[0]}{$id[1]} = $sequence;
}
close IN;

foreach my $id (keys %seq){
	my $seq = "";
	foreach my $index (sort keys %{$seq{$id}}){
		$seq = "$seq".$seq{$id}{$index};
	}
	print ">$id\n";
	print "$seq\n";
}
