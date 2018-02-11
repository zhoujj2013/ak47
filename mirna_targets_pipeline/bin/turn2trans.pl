#!/usr/bin/perl -w

use strict;

my $degs_f = shift;
my $refgene_f = shift;

my %g;
open IN,"$refgene_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my @trans = split /,/,$t[2];
	my @clean_trans;
	foreach my $tt (@trans){
		if($t[1] eq "protein coding gene"){
			if($tt =~ /NM/){
				push @clean_trans,$tt;
			}
		}else{
			push @clean_trans,$tt;
		}
	}
	$g{$t[0]} = \@clean_trans;
}
close IN;

open IN,"$degs_f" || die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	foreach my $trans_id (@{$g{$id}}){
		print "$trans_id\t$id\t";
		print join "\t",@t;
		print "\n";
	}
}
close IN;
