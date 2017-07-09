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

        This script create makefile for sequence annotation.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 <pep|dna> pretein.fa db.fa

USAGE
print "$usage";
exit(1);
};

my $type = shift;
my $fa = shift;
my $db = shift;

$db = abs_path($db);
$fa = abs_path($fa);

mkdir "./out" unless(-d "./out");
$c_outdir =  abs_path("./$out");

`ln -s $fa $c_outdir/protein.fa`;
$fa = "$c_outdir/protein.fa";

$mk .= "01cut_files.finished: $fa\n";
$mk .= "\tcd $c_outdir && perl $Bin/fastaDeal.pl --cutf 16 $fa && touch 01cut_files.finished\n";
$all .= "01cut_files.finished ";

$mk .= "02blast.finished: 01cut_files.finished\n";
$mk .= "\tcd $c_outdir && perl $Bin/create_shell.pl protein.fa.cut dna $db > blast.sh && perl $Bin/multi-process.pl -cpu 8 blast.sh > run.log 2>run.err && cat protein.fa.cut/*.outfmt6 > blast.outfmt6 && touch 02blast.finished\n";
$all .= "02blast.finished ";

$mk .= " ";
$mk .= "\t\n";
$all .= " ";

open OUT, ">$out/makefile";
print OUT $all, "\n";
print OUT $mk, "\n";
close OUT;
$all = "all: ";
$mk = "";

#########################

sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
