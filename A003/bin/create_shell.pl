#!/usr/bin/perl -w

use strict;

my $d = shift;
my $type = shift;
my $db = shift;

my @f = glob("$d/*");

foreach my $f (@f){
	if($type eq "pep"){
		print "blastp -query $f -db $db -out $f.outfmt6 -num_threads 1 -max_target_seqs 1 -outfmt 6 > $f.outfmt6.log 2>$f.outfmt6.err";
		print "\n";
	}elsif($type eq "dna"){
		print "blastx -query $f -db $db -out $f.outfmt6 -evalue 1e-20 -num_threads 1 -max_target_seqs 1 -outfmt 6 > $f.outfmt6.log 2>$f.outfmt6.err";
		print "\n";
	}
}
