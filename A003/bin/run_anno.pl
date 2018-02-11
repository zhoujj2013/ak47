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

        This script was designed to annotation protein or dna sequences by comparing to uniprot.
        Author: zhoujj2013\@gmail.com
        Usage: $0 query.fa dna/pep db.fa
        
USAGE
print "$usage";
exit(1);
};

my $q = shift;
my $type = shift;
my $db = shift;
my $prefix = shift;

my $db_name = basename($db, ".fasta");
my $db_dir = dirname($db);

my $bname = "";
if($q =~ /\.fa$/){
	$bname = basename($q, ".fa");
}elsif($q =~ /\.fasta$/){
	$bname = basename($q, ".fasta");
}else{
	print STDERR "Please input the query file in fasta format and named with .fa or .fasta.\n";
	exit(1);
}

#print "$bname\n";
`perl $Bin/fastaDeal.pl -cutf 16 $q`;
`perl $Bin/create_shell.pl ./$bname $type $db`;
`cat ./$bname/*.outfmt6 > ./blast.outfmt6`;
`perl $Bin/extract_uniprot.pl $db_dir/$db_name.dat ./blast.outfmt6 > $prefix.dbs.txt 2>$prefix.dbs.err`;
`perl $Bin/linkAnno.pl ./blast.outfmt6 $prefix.dbs.txt > $prefix.anno.txt 1>$prefix.anno.err`;

print STDERR "all tasks done!\n";



