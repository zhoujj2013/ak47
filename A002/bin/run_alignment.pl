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

        Run alignment for fastq files.
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

open OUT,">","$out/new_sample.lst" || die $!;
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;

	my $result_bam = "";
	
	my $phred = pop(@t);
	my $qual = "";
	if($phred eq "64"){
		$qual = "--phred64-quals";
	}
	
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
	
	mkdir "$out/$id" unless(-d "$out/$id");
	my $c_outdir =  "$out/$id";

	if($type eq "SE"){
		# tophat2
		if($conf{aligner} eq "tophat2"){
			$mk .= "$id.alignment.finished: $r1\n";
			$mk .= "\tcd $c_outdir && $conf{tophat2} $qual -G $conf{gtf} -p $conf{thread} -o ./ $conf{index} $r1 > $id.align.log 2> $id.align.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
			$result_bam = "$c_outdir/accepted_hits.bam";
		}
		# bowtie2	
		if($conf{aligner} eq "bowtie2"){
			$mk .= "$id.alignment.finished: $r1\n";
			$mk .= "\tcd $c_outdir && $conf{bowtie2} $qual -p $conf{thread} -x $conf{index} -U $r1 > $id.sam 2> $id.sam.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
			
			$mk .= "$id.sam2sortedbam.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && $conf{samtools} view -Sb ./$id.sam > $id.bam && $conf{samtools} sort -i $id.bam $id.sorted > sort.log 2>sort.err && cd - && touch $id.sam2sortedbam.finished\n";
			$all .= "$id.sam2sortedbam.finished ";

			$mk .= "$id.rmdup.finished: $id.sam2sortedbam.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}
		# bwa
		if($conf{aligner} eq "bwa"){
			$mk .= "$id.alignment.finished: $r1\n";
			$mk .= "\tcd $c_outdir && $conf{bwa} mem -t $conf{thread} -R '\@RG\\tID:foo\\tSM:$id\\tLB:$id.fq' -M $conf{index} $r1 > $id.sam 2> $id.sam.log && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
	
			$mk .= "$id.sam2sortedbam.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && $conf{samtools} view -Sb ./$id.sam > $id.raw.bam && $conf{samtools} fixmate $id.raw.bam $id.bam > fixmate.log 2>fixmate.err && $conf{samtools} sort -i $id.bam $id.sorted > sort.log 2>sort.err && cd - && touch $id.sam2sortedbam.finished\n";
			$all .= "$id.sam2sortedbam.finished ";

			$mk .= "$id.rmdup.finished: $id.sam2sortedbam.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}
		
		# star
		if($conf{aligner} eq "star"){
			$mk .= "$id.alignment.finished: $r1\n";
			$mk .= "\tcd $c_outdir && $conf{star} --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --runThreadN 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbGTFfile $conf{gtf} --alignIntronMax 1000000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --limitBAMsortRAM 4000000000  --genomeDir $conf{index} --readFilesIn $r1 --outFileNamePrefix ./ > $id.sam.log 2> $id.sam.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
	
			$mk .= "$id.rmdup.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}

	}elsif($type eq "PE"){
		# tophat2
		if($conf{aligner} eq "tophat2"){
			$mk .= "$id.alignment.finished: $r1 $r2\n";
			$mk .= "\tcd $c_outdir && $conf{tophat2} $qual -G $conf{gtf} -p $conf{thread} -o ./ $conf{index} $r1 $r2 > $id.align.log 2> $id.align.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
			$result_bam = "$c_outdir/accepted_hits.bam";
		}
		# bowtie2	
		if($conf{aligner} eq "bowtie2"){
			$mk .= "$id.alignment.finished: $r1 $r2\n";
			$mk .= "\tcd $c_outdir && $conf{bowtie2} $qual -p $conf{thread} -x $conf{index} -1 $r1 -2 $r2 > $id.sam 2> $id.sam.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
			
			$mk .= "$id.sam2sortedbam.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && $conf{samtools} view -Sb ./$id.sam > $id.bam && $conf{samtools} sort -i $id.bam $id.sorted > sort.log 2>sort.err && cd - && touch $id.sam2sortedbam.finished\n";
			$all .= "$id.sam2sortedbam.finished ";

			$mk .= "$id.rmdup.finished: $id.sam2sortedbam.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}
		# bwa
		if($conf{aligner} eq "bwa"){
			$mk .= "$id.alignment.finished: $r1 $r2\n";
			$mk .= "\tcd $c_outdir && $conf{bwa} mem -t $conf{thread} -R '\@RG\\tID:foo\\tSM:$id\\tLB:$id.fq' -M $conf{index} $r1 $r2 > $id.sam 2> $id.sam.log && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
	
			$mk .= "$id.sam2sortedbam.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && $conf{samtools} view -Sb ./$id.sam > $id.raw.bam && $conf{samtools} fixmate $id.raw.bam $id.bam > fixmate.log 2>fixmate.err && $conf{samtools} sort -i $id.bam $id.sorted > sort.log 2>sort.err && cd - && touch $id.sam2sortedbam.finished\n";
			$all .= "$id.sam2sortedbam.finished ";

			$mk .= "$id.rmdup.finished: $id.sam2sortedbam.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}
		
		# star
		if($conf{aligner} eq "star"){
			$mk .= "$id.alignment.finished: $r1\n";
			$mk .= "\tcd $c_outdir && $conf{star} --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --runThreadN 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbGTFfile $conf{gtf} --alignIntronMax 1000000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --limitBAMsortRAM 4000000000  --genomeDir $conf{index} --readFilesIn $r1 $r2 --outFileNamePrefix ./ > $id.sam.log 2> $id.sam.err && cd - && touch $id.alignment.finished\n";
			$all .= "$id.alignment.finished ";
	
			$mk .= "$id.rmdup.finished: $id.alignment.finished\n";
			$mk .= "\tcd $c_outdir && java -jar $conf{picard} MarkDuplicates I=./$id.sorted.bam O=./$id.sorted.rmdup.bam M=./$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ > ./$id.MarkDuplicates.log 2>./$id.MarkDuplicates.err && cd - && touch $id.rmdup.finished\n";
			$all .= "$id.rmdup.finished ";
			$result_bam = "$c_outdir/$id.sorted.rmdup.bam";
		}
	}
	
	#$mk .= "";
	#$mk .= "\t\n";
	#$all .= "$id.rmdup.finished ";
	
	print OUT "$id\t$result_bam\n";
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
