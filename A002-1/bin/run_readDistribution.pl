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

        Check reads distribution.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
$conf = abs_path($conf);

my %conf;
&load_conf($conf, \%conf);

my $out = abs_path($conf{OUTDIR});

my $all='all: ';
my $mk;

mkdir "$out" unless(-d "$out");

open OUT,">","$out/new_samples.lst" || die $!;
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	my $bam = $t[0];
	
	my $homer_path = "";
	
	mkdir "$out/$id" unless(-d "$out/$id");
	my $c_outdir =  "$out/$id";

	$mk .= "$id.bam2bed.finished: $bam\n";
	$mk .= "\tcd $c_outdir && $conf{bedtools} bamtobed -i $bam > ./$id.bed && cd - && touch $id.bam2bed.finished\n";
	$all .= "$id.bam2bed.finished ";
	
	$mk .= "$id.makeDir.finished: $id.bam2bed.finished\n";
	$mk .= "\tcd $c_outdir && $conf{homer}/makeTagDirectory $id.homer $c_outdir/$id.bed -genome $conf{spe} -checkGC -format bed >$id.makeDir.log 2> $id.makeDir.err && cd - && touch $id.makeDir.finished\n";
	$all .= "$id.makeDir.finished ";
	
	$mk .= "$id.dist.finished: $id.makeDir.finished\n";
	$mk .= "\tcd $c_outdir && perl $conf{homer}/annotatePeaks.pl $c_outdir/$id.bed $conf{spe} > /dev/null 2>$id.dist.log && cd - && touch $id.dis.finished\n";
	$all .= "$id.dist.finished ";
	
	$homer_path = "$c_outdir/$id.homer";
	print OUT "$id\t$homer_path\n";
}
close IN;
close OUT;

#### write you things ###
make_path abs_path($conf{OUTDIR});
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
