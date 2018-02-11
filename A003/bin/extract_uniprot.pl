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
	my $id = $1 if($con =~ /ID   (\S+)/);
	my $acc = $1 if($con =~ /\nAC   ([^;]+);/);
	next if($acc eq "");
	next unless(exists $sw{$acc});
	
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
		if($l =~ /DR   InterPro; ([^;]+); (.*);/){
			push @ipr,"$1~$2";
		}
		if($l =~ /DR   Pfam; ([^;]+); (.*);/){
			push @pfam,"$1~$2";
		}
	}
	
	my %g;
	if(scalar(@go) > 0){
		foreach my $g (@go){
			my ($aspect,$d) = ($1,$2) if($g->[1] =~ /^([P|F|C]):([^;]+);/);
			#print "$aspect\t$d\n";
			if($aspect eq "F"){
				push @{$g{MF}},"$g->[0]~$d";
			}elsif($aspect eq "P"){
				push @{$g{BP}},"$g->[0]~$d";
			}elsif($aspect eq "C"){
				push @{$g{CC}},"$g->[0]~$d";
			}
		}
	}
	
	#print Dumper(\@go);
	#print Dumper(\@ipr);
	#print Dumper(\@pfam);
	if($refseq eq ""){
		$refseq = "NA";
	}
	
	if($kegg eq ""){
		$kegg = "NA";
	}else{
		my $kegg_id = $kegg;
		sleep(0.5);
		my $kegg_con = `wget -q http://rest.kegg.jp/get/$kegg`;
		my @kegg_con = split /\n/,$kegg_con;
		my $flag = 0;
		my @p;
		my $k_str = "";
		open KEGG,"./$kegg" || die $!;
		$/ = "\n";
		while(<KEGG>){
			chomp;
			my $kl = $_;
			#print STDERR "$kl\n";
			if($kl =~ /^ORTHOLOGY\s+(\S+)\s+(.*)/){
				my $k_id = $1;
				my $k_desc = $2;
				$k_str = "$1~$2";
			}
			
			if($kl =~ /^PATHWAY\s+(\S+)\s+(.*)/){
				my $p_id = $1;
				my $p_desc = $2;
				$p_id = $1 if($p_id =~ /(\d+)$/);
				$p_id = "map$p_id";
				push @p,"$p_id~$p_desc";
				$flag = 1;
				next;
			}

			if($kl =~ /^BRITE/){
				$flag = 0;
			}
	
			if($kl =~ /\s+(\S+)\s+(.*)/ && $flag == 1){
				my $p_id = $1;
				my $p_desc = $2;
				$p_id = $1 if($p_id =~ /(\d+)$/);
				$p_id = "map$p_id";
				push @p,"$p_id~$p_desc";
			}
		}
		close KEGG;
		$/ = "//";
		
		if($k_str eq ""){
			$k_str = "NA";
		}
		
		my $str;
		if(scalar(@p) < 1){
			$str = "NA";
		}else{
			$str = join "; ", @p;
		}
		$kegg = "$kegg\t$k_str\t$str";
		`rm $kegg_id`;
	}
	
	print "$id\t$acc\t$refseq\t$kegg\t";

	# output GO annotation
	my @go_type = ("BP","MF","CC");
	foreach my $go_t (@go_type){
		#print "$go_t\n";
		if(exists $g{$go_t}){
			print join "; ",@{$g{$go_t}};
			print "\t";
		}else{
			print "NA\t";
		}
	}
	
	# output IPR annotation
	if(scalar(@ipr) > 0){
		print join "; ",@ipr;
		print "\t";
	}else{
		print "NA\t";
	}

	# output pfam
	if(scalar(@pfam) > 0){
		print join "; ",@pfam;
	}else{
		print "NA";
	}
	
	print "\n";
}
close IN;
$/ = "\n";


