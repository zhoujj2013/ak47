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

        Create tables for report.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jul. 8 16:06:43 HKT 2017
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

mkdir "./result" unless(-d "./result");

#### Table 1
my %qc;
my @id;
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	push @id,$id;
	my $type = "";
	if(scalar(@t) > 1){
		$type = "PE";
	}else{
		$type = "SE";
	}

	#print "$type\n";
	my @arr;
	#print "$out/$id/Trimmomatic.err\n";
	open TRIM,"$out/$id/Trimmomatic.err" || die $!;
	while(<TRIM>){
		chomp;
		#print $_,"\n";
		next unless(/^Input Read/);
		#print $_,"\n";
		#Input Read Pairs: 10000 Both Surviving: 5527 (55.27%) Forward Only Surviving: 2440 (24.40%) Reverse Only Surviving: 1657 (16.57%) Dropped: 376 (3.76%)
		#Input Reads: 10000 Surviving: 6481 (64.81%) Dropped: 3519 (35.19%)
		if($type eq "PE"){
			my $all_reads = $1*2 if(/Input Read Pairs: (\d+)/);
			my $surviving = $1 if(/Both Surviving: (\d+)/);
			my $forward_only = $1 if(/Forward Only Surviving: (\d+)/);
			my $reverse_only = $1 if(/Reverse Only Surviving: (\d+)/);
			my $dropped = $1 if(/Dropped: (\d+)/);
			my $clean_reads = $surviving*2 + $forward_only + $reverse_only;
			my $filtered_reads = $dropped*2;
			my $clean_reads_perc = sprintf("%.2f",$clean_reads*100/$all_reads);
			@arr = ($all_reads, $surviving, $clean_reads, $filtered_reads, $clean_reads_perc);
		}elsif($type eq "SE"){
			my $all_reads = $1 if(/Input Reads: (\d+)/);
			my $surviving = $1 if(/Surviving: (\d+)/);
			my $dropped = $1 if(/Dropped: (\d+)/);
			my $clean_reads = $surviving;
			my $filtered_reads = $dropped;
			my $clean_reads_perc = sprintf("%.2f",$clean_reads*100/$all_reads);
			@arr = ($all_reads,"NA",$clean_reads,$filtered_reads, $clean_reads_perc);
		}
	}
	close TRIM;
	#print @arr,"\n";
	my $gc = 0;
	if($type eq "PE"){
		my $gc1 = `grep "^%GC" out/$id/$id.trimmed.rd.R1.fq_fastqc/fastqc_data.txt | cut -f 2`;
		my $gc2 = `grep "^%GC" out/$id/$id.trimmed.rd.R2.fq_fastqc/fastqc_data.txt | cut -f 2`;	
		chomp($gc1);
		chomp($gc2);
		$gc = ($gc1 + $gc2)/2;
	}elsif($type eq "SE"){
		$gc = `grep "^%GC" out/$id/$id.trimmed.rd.R1.fq_fastqc/fastqc_data.txt | cut -f 2`;
		chomp($gc);
	}
	push @arr,$gc;
	$qc{$id} = \@arr;
}
close IN;

open OUT,">","./result/table1.xls" || die $!;
print OUT "Table 1\n";
print OUT "Sample\tRaw Reads\tClean Paired reads\tClean reads\tFiltered reads\tClean reads(%)\tGC(%)\n";
print "##########Table 1##########\n";
print "|Sample |Raw Reads|Clean Paired reads|Clean reads|Filtered reads|Clean reads(%)|GC(%)|\n";
print "| -------- | :-----:  | :----:  | :----: | :----: |:----: |:----:|\n";
foreach my $id (@id){
	print OUT "$id\t";
	print OUT join "\t",@{$qc{$id}};
	print OUT "\n";
	print "|$id|";
	print join "|",@{$qc{$id}};
	print "|\n";
}
print "#########################\n\n\n\n";
close OUT;

#### Table2
my @dup;
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	my $too_short = `grep "^Short read :" $out/$id/$id.trimmed.rd.log | awk '{print \$4}'`;
	my $N_in_seq = `grep "^N in seq" $out/$id/$id.trimmed.rd.log | awk '{print \$5}'`;
	my $Valid_read = `grep "^Valid read" $out/$id/$id.trimmed.rd.log | awk '{print \$4}'`;
	my $repeat = `grep "^Repeat CC" $out/$id/$id.trimmed.rd.log | awk '{print \$4}'`;
	my $Uniq_read = `grep "^Uniq read" $out/$id/$id.trimmed.rd.log | awk '{print \$4}'`;
	my $Duplicate = `grep "^Duplicate" $out/$id/$id.trimmed.rd.log | awk '{print \$3}'`;

	chomp($too_short);
	chomp($N_in_seq);
	chomp($Valid_read);
	chomp($repeat);
	chomp($Uniq_read);
	chomp($Duplicate);
	
	my $clean_reads = $too_short + $N_in_seq + $Valid_read + $repeat + $Uniq_read + $Duplicate;
	my $dup_rate = sprintf("%.2f",$Duplicate*100/$clean_reads);
	push @dup, [$id,$clean_reads,$Duplicate,$Uniq_read,$dup_rate];
}
close IN;

open OUT,">","./result/table2.xls" || die $!;
print OUT "Table 2\n";
print OUT "Sample\tClean reads\tDup. reads\tUniq. reads\tDup. rate(%)\n";
print "####### Table 2 ###########\n";
print "|Sample |Clean reads|Dup. reads|Uniq. reads|Dup. rate(%)|\n";
print "| --------   | :-----:  | :----:  | :----: |:----: |\n";
foreach my $d (@dup){
	print OUT join "\t",@$d;
	print OUT "\n";
	print "|";
	print join "|",@$d;
	print "|\n";
	
}
print "###########################\n\n\n\n";
close OUT;

# create data set
open IN,"$conf{SAMPLE}" || die $!;
while(<IN>){
	chomp;
    my @t = split /\t/;
    my $id = shift @t;
    push @id,$id;
    my $type = "";
    if(scalar(@t) > 1){
        $type = "PE";
    }else{
        $type = "SE";
    }
	
	mkdir "./result/$id" unless(-d "./result/$id");

	if($type eq "PE"){
		`cp $out/$id/$id.trimmed.R1.fq ./result/$id/`;
		`cp $out/$id/$id.trimmed.R2.fq ./result/$id/`;
		`cp $out/$id/$id.trimmed.R1.unpaired.fq ./result/$id/`;
		`cp $out/$id/$id.trimmed.R2.unpaired.fq ./result/$id/`;
		
		`cp $out/$id/$id.trimmed.rd.R1.fq ./result/$id/`;
		`cp $out/$id/$id.trimmed.rd.R2.fq ./result/$id/`;

		`cp -r $out/$id/$id.trimmed.rd.R1.fq_fastqc ./result/$id/`;
		`cp -r $out/$id/$id.trimmed.rd.R2.fq_fastqc ./result/$id/`;
	}elsif($type eq "SE"){
		`cp $out/$id/$id.trimmed.R1.fq ./result/$id/`;
		
		`cp $out/$id/$id.trimmed.rd.R1.fq ./result/$id/`;

		`cp -r $out/$id/$id.trimmed.rd.R1.fq_fastqc ./result/$id/`;
	}
}
close OUT;

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
        #warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
