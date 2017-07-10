#!/usr/bin/perl 

use strict;
use Data::Dumper;

my $f = shift;
my $blast = shift;

#print $f,"\n";

# read in blast file
my %b;
my %sw;
open IN,"$blast" || die $!;
while(<IN>){
	chomp;
	#print "$_";
	my @t = split /\t/;
	$b{$t[0]} = \@t;
	my @id = split /\|/,$t[1];
	my $acc = $id[1];
	push @{$sw{$acc}}, \@t;;
}
close IN;

#print Dumper(\%b);

# read in swissprot file
my %acc;
open IN,"$f" || die $!;
$/ = "//";
while(<IN>){
	chomp;
	my $con = $_;
	#print $con,"\n";
	my @l = split /\n/,$con;
	my $acc = $1 if($con =~ /\nAC   ([^;]+);/);
	next if($acc eq "");
	#next unless(exists $sw{$acc});
	
	my $desc = "";
	my @go;
	my $kegg;
	my @ipr;
	my $geneid = "";
	my $refseq = "";
	my @pfam;
	
	foreach my $l (@l){
		if($l =~ /^DR   RefSeq; ([^;]+);/){
			$refseq = $1;
		}
		if($l =~ /^DR   GeneID; ([^;]+);/){
			$geneid = $1;
		}
		if($l =~ /DR   KEGG; ([^;]+);/){
			$kegg = $1;
		}
		if($l =~ /DR   GO; ([^;]+); (.*)/){
			push @go,[$1, $2];
		}
		if($l =~ /DR   InterPro; ([^;]+); (.*)/){
			push @ipr,[$1, $2];
		}
		if($l =~ /DR   Pfam; ([^;]+); (.*)/){
			push @pfam,[$1, $2];
		}
	}
	#print Dumper(\@go);
	#print Dumper(\@ipr);
	#print Dumper(\@pfam);
	print $acc,"\t$refseq\t$kegg\t@go\n";
}
close IN;
$/ = "\n";
