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

        This script was designed to link annotation to novel sequences.
        Author: zhoujj2013\@gmail.com
        Usage: $0 blast.outfmt6 dbs.txt > result.txt
        
USAGE
print "$usage";
exit(1);
};

my $b = shift;
my $dbs = shift;

my %anno;
my $col_num;
open IN,"$dbs" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[0];
	$anno{$id} = \@t;
	$col_num = scalar(@t);
}
close IN;

my @na;
for(my $i = 0; $i < $col_num; $i++){
	push @na,"NA";
}

open IN,"$b" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $novel_id = $t[0];
	my $sw_id = $t[1];
	my @sw_id = split /\|/,$sw_id;
	if(exists $anno{$sw_id[2]}){
		print "$novel_id\t";
		print join "\t",@{$anno{$sw_id[2]}};
		print "\n";
	}else{
		print "$novel_id\t";
		print join "\t",@na;
		print "\n";
	}
}
close IN;

