#!/usr/bin/perl -w -s

#use lib "/Users/d.grun/data/bin";
# use lib "/home/d.grun/bin";
use lib "/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation";
use tools;


if (scalar @ARGV == 1){
    die "usage: -in=INPUT.gtf -ref=GENOME.fa -f=ID field name\n" if $ARGV[0] eq "help";
}
$f = "transcript_id" if !$f;

%genome = ();
fasta2hash(\%genome,$ref);
$j = 0;
$flag  = 0;
%first = ();
open(IN,"<",$in);
while(<IN>){
  chomp;
  @F = split(/\t/);
  next if $F[2] =~ /transcript/;
  if ( !exists($genome{$F[0]}) ){
    $ns{$F[0]} = 1;
    next;
  }
  @G = split(/\s/,$F[$#F]);
  $id = get_id($_,$f);
#  $id = $G[3];
#  $id =~ s/[\";]//g;
  if (!exists($first{$id})){
    if ($flag){
      $seq = get_sequence(\%exons,\%genome,$region);
      $seq = revcompl($seq) if $str eq "-";
      print ">".$ID."\n".$seq."\n";
    }
    $j ++;
    print STDERR $j."\r";
    $flag       = 1;
    $first{$id} = 1;
    $ID         = $id;
    $region     = $F[0];
    $str        = $F[6];
    %exons      = ();
  }
  $exons{$F[3]} = $F[4];
}
$seq = get_sequence(\%exons,\%genome,$region);
$seq = revcompl($seq) if $str eq "-";
print ">".$id."\n".$seq."\n" if exists($genome{$region});
close(IN);

foreach (sort keys %ns ){
  print STDERR $_."\n";
}

sub get_sequence {
  $qhash = shift;
  $rhash = shift;
  $key   = shift;
  $flag  = 1;
  foreach $k (sort {$a<=>$b} keys %$qhash){
    if ($flag){
      $seq = substr($$rhash{$key}, $k-1, $$qhash{$k} - $k + 1);
    }else{
      $seq = $seq.substr($$rhash{$key}, $k-1, $$qhash{$k} - $k + 1);
    }
    $flag = 0;
  }
  return $seq;
}

sub get_id {
  $x = shift;
  $f = shift;
  chomp($x);
  %h = ();
  @H = split(/\t/,$_);
  $H[$#H] =~ s/[\";]//g;
  @G = split(/\s/,$H[$#H]);

  $i = 0;
  while ( $i < $#G ){
    $key = $G[$i];
    $i++;
    $value = $G[$i];
    $h{$key} = $value;
  }
  return($h{$f});
}
