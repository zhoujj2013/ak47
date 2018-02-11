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

my $target_f = shift;
my $f=shift;

my %target;
open IN,"$target_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]";
	$target{$id} = 1;
}
close IN;

open IN,"$f" || die $!;
$/ = ">"; <IN>; $/ = "\n";
while(<IN>){
	chomp;
	my $header = $_;
	my ($mir_id, $trans_id) = ($1,$2) if(/(\S+)\s+(\S+)/);
	$/ = ">";
	my $con = <IN>;
	chomp($con);
	$/ = "\n";
	my $id = "$mir_id#$trans_id";
	if(exists $target{$id}){
		print ">$header\n";
		print "$con";
	}
}
close IN;
