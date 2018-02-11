#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;

&usage if @ARGV<1;

#open IN,"" ||die "Can't open the file:$\n";
#open OUT,"" ||die "Can't open the file:$\n";

sub usage {
        my $usage = << "USAGE";

        Description of this script. 5/21/2016
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xxx.cfg

USAGE
print "$usage";
exit(1);
};

my ($deg_f, $expr_f, $fc_cutoff) = @ARGV;

open OUT1,">","./up.fc.txt" || die $!;
open OUT2,">","./down.fc.txt" || die $!;
my %deg;
open IN,"$deg_f" || die $!;
my $header = <IN>;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $fc = $t[5];
	$deg{$t[0]} = \@t if($fc >= $fc_cutoff || $fc <= -$fc_cutoff);
	print OUT1 "$_\n" if($fc >= $fc_cutoff);
	print OUT2 "$_\n" if($fc <= -$fc_cutoff);
	#print $fc,"\n"  if($fc >= 1 || $fc <= -1);
}
close IN;
close OUT1;
close OUT2;

open IN,"$expr_f" || die $!;
my $h = <IN>;
print "$h";
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[0];
	if(exists $deg{$id}){
		print join "\t",@t;
		print "\n";
	}
}
close IN;
