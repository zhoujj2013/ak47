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

        Clean fastq dataset.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 config.cfg
        default(rnaseq parameter, rm duplications, minlen >= 36)

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
if($conf eq "default" || $conf eq "0"){
	$conf = abs_path("$Bin/config.txt");
}else{
	$conf = abs_path($conf);
}
my %conf;
&load_conf($conf, \%conf);

my $out = abs_path($conf{OUTDIR});

my $all='all: ';
my $mk;



mkdir "$out" unless(-d "$out");

open OUT,">","$out/new_sample.lst" || die $!;
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;

	my $cleaned_r1 = "";
	my $cleaned_r2 = "";
	
	my $type = "";
	my $r1 = "";
	my $r2 = "";
	if(scalar(@t) == 1){
		$type = "SE";
		$r1 = abs_path($t[0]);
	}else{
		$type = "PE";
		$r1 = abs_path($t[0]);
		$r2 = abs_path($t[1]);
	}
	my $phred = `perl $Bin/phredDetector.pl $r1`;
	
	mkdir "$out/$id" unless(-d "$out/$id");
	my $c_outdir =  "$out/$id";

	if($type eq "SE"){
		$mk .= "$id.fastqc.finished: $r1\n";
		$mk .= "\tcd $c_outdir && $conf{fastqc} $r1 -o ./ > $id.fastqc.1st.log 2> $id.fastqc.1st.err && cd -&& touch $id.fastqc.finished\n";
		$all .= "$id.fastqc.finished ";
		
		$mk .= "$id.trimming.finished: $r1\n";
		$mk .= "\tcd $c_outdir && java -jar $conf{Trimmomatic} SE -phred$phred $r1 $id.trimmed.R1.fq ILLUMINACLIP:$conf{ILLUMINACLIP}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > ./Trimmomatic.log 2>./Trimmomatic.err && cd - && touch $id.trimming.finished\n";
		$all .= "$id.trimming.finished ";
		
		if($conf{LIB} eq "rnaseq"){
			$mk .= "$id.rmdup.finished: $id.trimming.finished\n";
			$mk .= "\tcd $c_outdir && $conf{rmdup} -a ./$id.trimmed.R1.fq -b ./$id.trimmed.R1.fq -o $id.trimmed.rd > $id.rmdup.log 2>$id.rmdup.err && rm $id.trimmed.rd.R2.fq && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";

			$mk .= "$id.fastqc.2nd.finished: $id.rmdup.finished\n";
			$mk .= "\tcd $c_outdir && $conf{fastqc} $id.trimmed.rd.R1.fq -o ./ > $id.fastqc.2nd.log 2> $id.fastqc.2nd.err && cd -&& touch $id.fastqc.2nd.finished\n";
			$all .= "$id.fastqc.2nd.finished ";
			$cleaned_r1 = "$out/$id/$id.trimmed.rd.R1.fq";
		}elsif($conf{LIB} eq "dnaseq"){
			$mk .= "$id.rmdup.finished: $id.trimming.finished\n";
            $mk .= "\tcd $c_outdir && echo \"Don't run remove duplications!\" > $id.rmdup.log 2>$id.rmdup.err && cd - && touch $id.rmdup.finished\n";
            $all .= "$id.rmdup.finished ";
	
			$mk .= "$id.fastqc.2nd.finished: $id.rmdup.finished\n";
			$mk .= "\tcd $c_outdir && $conf{fastqc} $id.trimmed.rd.R1.fq -o ./ > $id.fastqc.2nd.log 2> $id.fastqc.2nd.err && cd -&& touch $id.fastqc.2nd.finished\n";
			$all .= "$id.fastqc.2nd.finished ";
			$cleaned_r1 = "$out/$id/$id.trimmed.rd.R1.fq";
		}
	}elsif($type eq "PE"){
		#fastqc
		$mk .= "$id.fastqc.finished: $r1 $r2\n";
		$mk .= "\tcd $c_outdir && $conf{fastqc} $r1 $r2 -o ./ > $id.fastqc.1st.log 2> $id.fastqc.1st.err && cd -&& touch $id.fastqc.finished\n";
		$all .= "$id.fastqc.finished ";
		#trimming
		$mk .= "$id.trimming.finished: $r1 $r2\n";
		$mk .= "\tcd $c_outdir && java -jar $conf{Trimmomatic} PE -phred$phred $r1 $r2 $id.trimmed.R1.fq $id.trimmed.R1.unpaired.fq $id.trimmed.R2.fq $id.trimmed.R2.unpaired.fq ILLUMINACLIP:$conf{ILLUMINACLIP}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > ./Trimmomatic.log 2>./Trimmomatic.err && cd - && touch $id.trimming.finished\n";
		$all .= "$id.trimming.finished ";
		
		if($conf{LIB} eq "rnaseq"){
			$mk .= "$id.rmdup.finished: $id.trimming.finished\n";
			$mk .= "\tcd $c_outdir && $conf{rmdup} -a $id.trimmed.R1.fq -b $id.trimmed.R2.fq -o $id.trimmed.rd > $id.rmdup.log 2>$id.rmdup.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";

			$mk .= "$id.fastqc.2nd.finished: $id.rmdup.finished\n";
			$mk .= "\tcd $c_outdir && $conf{fastqc} $id.trimmed.rd.R1.fq $id.trimmed.rd.R2.fq -o ./ > $id.fastqc.2nd.log 2> $id.fastqc.2nd.err && cd -&& touch $id.fastqc.2nd.finished\n";
			$all .= "$id.fastqc.2nd.finished ";
			$cleaned_r1 = "$out/$id/$id.trimmed.rd.R1.fq";
			$cleaned_r2 = "$out/$id/$id.trimmed.rd.R2.fq";
		}elsif($conf{LIB} eq "dnaseq"){
			$mk .= "$id.rmdup.finished: $id.trimming.finished\n";
            $mk .= "\tcd $c_outdir && echo \"Don't run remove duplications!\" > $id.rmdup.log 2>$id.rmdup.err && cd - && touch $id.rmdup.finished\n";
            $all .= "$id.rmdup.finished ";
	
			$mk .= "$id.fastqc.2nd.finished: $id.rmdup.finished\n";
			$mk .= "\tcd $c_outdir && $conf{fastqc} $id.trimmed.R1.fq $id.trimmed.R2.fq -o ./ > $id.fastqc.2nd.log 2> $id.fastqc.2nd.err && cd -&& touch $id.fastqc.2nd.finished\n";
			$all .= "$id.fastqc.2nd.finished ";
			$cleaned_r1 = "$out/$id/$id.trimmed.R1.fq";
            $cleaned_r2 = "$out/$id/$id.trimmed.R2.fq";
		}
	}
	print OUT "$id\t$cleaned_r1";
	if($cleaned_r2 eq ""){
		print OUT "\t$phred\n";
	}elsif($cleaned_r2 ne ""){
		print OUT "\t$cleaned_r2\t$phred\n";
	}
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
