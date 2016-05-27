#!/usr/bin/perl -w -s

use lib "/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation";  	# use lib '/hpc/hub_oudenaarden/bin/';
use tools;


if (scalar @ARGV == 1){
    die "usage: -in=input.fa -n=number of As (default: 10) -out=output.fa\n" if $ARGV[0] eq "help";
}

$n = 10 if !$n;
%seqs = ();
fasta2hash(\%seqs,$in);
$m = "(A|a){".$n.",}";

open(OUT,">",$out);
foreach $k (keys %seqs){
  $s = $seqs{$k};
  while ( $s =~ /$m/g ){
    $x = $&;
    $x =~ s/[Aa]/N/g;
    $s = substr($s,0,pos($s) - length($x)).$x.substr($s,pos($s),length($s) - pos($s));
  }
  print OUT ">".$k."\n".$s."\n";
}
close(OUT);

