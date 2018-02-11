#!/usr/bin/perl -w

use strict;

my $f = shift;

open IN,"$f" || die $!;
<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;<IN>;
$/ = "Complete";
while(<IN>){
	chomp;
	my $con = $_;
	chomp($con);
	#print $con;
	my @l = split /\n/,$con;
	my @pos;
	my $hit_flag = 0;
	foreach my $l (@l){
		if($l =~ /No Hits Found above Threshold/){
			$hit_flag = 1;
		}
		if($l =~ /^\s+/){
			next if($l =~ /^\s+Forward:/);
			push @pos,$l;
		}elsif($l =~ /^>\w+/){
			push @pos,$l;
		}
	}
	next if($hit_flag == 1);
	
	my @l_tmp;
	foreach my $p (@pos){
		if($p =~ />/){
			print "$p\n";
			print join "\n",@l_tmp;
			print "\n";
			@l_tmp = "";
		}else{
			push @l_tmp,$p;
		}
		#print "$p\n";
	}
}
close IN;
