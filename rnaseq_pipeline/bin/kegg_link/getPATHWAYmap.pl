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
        Author: zhoujj2013\@gmail.com 
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

my $up_f=shift;
my $down_f=shift;
my $geneid_f=shift;

my %geneid;
open IN,"$geneid_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$geneid{$t[0]} = $t[1];
}
close IN;

my %kegg_c;
my %up;
open IN,"$up_f" || die $!;
while(<IN>){
	chomp;
	next unless(/^KEGG_PATHWAY/);
	my @t = split /\t/;
	my $kegg = $t[11];
	my $gene = $t[12];
	#print "$kegg\n";
	my @gene = split /, /,$gene;
	my @geneid;
	foreach my $g (@gene){
		next unless($geneid{$g});
		push @geneid,$geneid{$g};
	}
	my $gene_str = join "+red%0d%0a",@geneid;
	my @kegg = split /:/,$kegg;
	# http://www.kegg.jp/kegg-bin/show_pathway?map=mmu05144&multi_query=16176+red%0d%0a21926+green
	$kegg_c{$kegg[0]} = $kegg;
	#my $out = "$t[11]\t";
	#$out .= "http://www.kegg.jp/kegg-bin/show_pathway?map=";
	#$out .= "$kegg[0]";
	#my $out = "&multi_query=";
	my $out = "$gene_str+red";
	$up{$kegg[0]} = $out;
	#print "\n";
}
close IN;

my %down;
open IN,"$down_f" || die $!;
while(<IN>){
	chomp;
	next unless(/^KEGG_PATHWAY/);
	my @t = split /\t/;
	my $kegg = $t[11];
	my $gene = $t[12];
	#print "$kegg\n";
	my @gene = split /,\s+/,$gene;
	my @geneid;
	foreach my $g (@gene){
		next unless($geneid{$g});
		push @geneid,$geneid{$g};
	}
	my $gene_str = join "+cyan%0d%0a",@geneid;
	my @kegg = split /:/,$kegg;
	# http://www.kegg.jp/kegg-bin/show_pathway?map=mmu05144&multi_query=16176+red%0d%0a21926+green
	$kegg_c{$kegg[0]} = $kegg;
	#my $out = "$t[11]\t";
	#$out .= "http://www.kegg.jp/kegg-bin/show_pathway?map=";
	#$out .= "$kegg[0]";
	#my $out = "&multi_query=";
	my $out = "$gene_str+cyan";
	$down{$kegg[0]} = $out;
	#print "\n";
}
close IN;

foreach my $id (keys %kegg_c){
	my @id = split /:/,$kegg_c{$id};
	my $query = "";
	if(!(exists $down{$id} && exists $up{$id})){
		if(exists $down{$id}){
			$query = $down{$id};
		}
		if(exists $up{$id}){
			$query = $up{$id};
		}
	}elsif(exists $down{$id} && exists $up{$id}){
		$query = $down{$id}."%0d%0a".$up{$id};
	}
	my $link = "http://www.kegg.jp/kegg-bin/show_pathway?map=";
	$link .= "$id";
	$link .= "&multi_query=";
	$link .= "$query";
	print "$kegg_c{$id}\t$link\n";
}
