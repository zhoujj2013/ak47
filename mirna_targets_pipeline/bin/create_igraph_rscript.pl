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

        This script create makefile for LncFunNet analysis.
        Version: v1.0
        Author: zhoujj2013\@gmail.com
        Last modified: Wed Jun 14 16:06:43 HKT 2017
        Usage: $0 config.cfg
        
        NOTE: please check config.cfg format in ./demo directory.

USAGE
print "$usage";
exit(1);
};

my $interaction_f = shift;
my $top_num = shift;
my $prefix = shift;

my %int;
open IN,"$interaction_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	print STDERR "$t[0]\n" unless(defined $t[3]);
	$int{$t[0]}{$t[3]} = 1;
	#if($i < $top_num){
	#	push @mir,$t[0];
	#	push @gene,$t[3];
	#}
}
close IN;

my $i = 0;
my @mir;
my @gene;
foreach my $m (keys %int){
	my $m2 = "\"$m\"";
	if($i < $top_num){
		foreach my $g (keys %{$int{$m}}){
			$g = "\"$g\"";
			push @mir,$m2;
			push @gene,$g;
		}
	}else{
		last;
	}
	$i++;	
}

my $mir = join ",",@mir;
#print "$mir\n";
my $gene = join ",",@gene;
my %node;
foreach my $m (@mir){
	$node{$m} = '"orange"';
}
foreach my $g (@gene){
	$node{$g} = '"lightgreen"';
}

my @node_name;
my @node_color;
foreach my $n (keys %node){
	push @node_name, "$n";
	push @node_color, "$node{$n}";
}

my $node_name_str = join ",",@node_name;
my $node_color_str = join ",",@node_color;

print "library(\"igraph\")\n";
print "nodes = data.frame(name=c($node_name_str), color=c($node_color_str))\n";
print "relations <- data.frame(from=c($mir), to=c($gene))\n";
print "g <- graph.data.frame(relations, vertices=nodes)\n";
print "V(g)\$color = vertex_attr(g, \"color\")\n";
print "pdf(file=\"$prefix.pdf\")\n";
print "plot(g, edge.arrow.size=0.5, layout=layout_with_fr, vertex.size=10)\n";
print "dev.off()\n";


