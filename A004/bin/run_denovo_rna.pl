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

        This script designed to run trinity pipeline.
        Author: zhoujj2013\@gmail.com
        Usage: $0 r1.fq r2.fq prefix
        
USAGE
print "$usage";
exit(1);
};

my $prefix = pop(@ARGV);
my @fq = @ARGV;

my $Trinity = "/home/zhoujj/software/trinity/";

### assembly ##########
my ($r1, $r2);
if(scalar(@fq) == 1){
	$r1 = abs_path($ARGV[0]);
	`$Trinity/Trinity --seqType fq --max_memory 15G --single $r1 --CPU 4 > trinity.log 2> trinity.err`;
}elsif(scalar(@fq) == 2){
	$r1 = $ARGV[0];
	$r2 = $ARGV[1];
	`$Trinity/Trinity --seqType fq --max_memory 15G --left $r1 --right $r2 --CPU 4 > trinity.log 2> trinity.err`;
}

`perl $Trinity/util/TrinityStats.pl ./trinity_out_dir/Trinity.fasta > $prefix.contig.Nx.stat`;

## recheck the reads
if(scalar(@fq) == 1){
    $r1 = abs_path($ARGV[0]);
	`perl $Trinity/util/align_and_estimate_abundance.pl --transcripts ./trinity_out_dir/Trinity.fasta --seqType fq --single $r1 --est_method RSEM --aln_method bowtie2 --trinity_mode --output_dir rsem_outdir --prep_reference > rsem.log 2> rsem.err`;
}elsif(scalar(@fq) == 2){
	$r1 = abs_path($ARGV[0]);
    $r2 = abs_path($ARGV[1]);
	`perl $Trinity/util/align_and_estimate_abundance.pl --transcripts ./trinity_out_dir/Trinity.fasta --seqType fq --left $r1 --right $r2 --est_method RSEM --aln_method bowtie2 --trinity_mode --output_dir rsem_outdir --prep_reference > rsem.log 2> rsem.err`;
}

######### Annotate protein sequences ##############
my $TransDecoder = "/home/zhoujj/software/TransDecoder-3.0.1";
my $db = "/disk4/database/uniprot/uniprot_sprot.fasta";

`$TransDecoder/TransDecoder.LongOrfs -t ./trinity_out_dir/Trinity.fasta > TransDecoder.LongOrfs.log 2>TransDecoder.LongOrfs.err`;

`perl $Bin/fastaDeal.pl -cutf 16 Trinity.fasta.transdecoder_dir/longest_orfs.pep`;
open OUT,">","blast.sh" || die $!;
foreach my $f (glob("./longest_orfs.pep.cut/*")){
	print OUT "blastp -query $f  -db $db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > $f.outfmt6\n";
}
close OUT;
`perl $Bin/multi-process.pl -cpu 8 blast.sh > blast.log 2>blast.err`;
`cat longest_orfs.pep.cut/*.outfmt6 > blast.outfmt6`;

my $hmmscan = "/home/zhoujj/software/hmmer/binaries/hmmscan";
my $pfam = "/disk4/database/pfam/Pfam-A.hmm";
#hmmscan --cpu 8 --domtblout pfam.domtblout /path/to/Pfam-A.hmm transdecoder_dir/longest_orfs.pep
`$hmmscan --cpu 8 --domtblout pfam.domtblout $pfam ./Trinity.fasta.transdecoder_dir/longest_orfs.pep > pfam.log 2>pfam.err`;

`$TransDecoder/TransDecoder.Predict -t ./trinity_out_dir/Trinity.fasta --retain_pfam_hits ./pfam.domtblout --retain_blastp_hits ./blast.outfmt6 > TransDecoder.Predict.log 2>TransDecoder.Predict.err`;


