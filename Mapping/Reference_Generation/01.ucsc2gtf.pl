#!/usr/bin/perl -w -s

#use lib "/Users/d.grun/data/bin";
# use lib "/home/d.grun/bin";
use lib "/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation";
use tools;

if (scalar @ARGV == 1)
{
    die "usage: -in=INPUT.ucsc_format -out=OUTPUT.gtf -m=gene2isoforms.tsv -utr=0,3,5\n" if ($ARGV[0] eq "help");
}
$utr = 0 if !$utr;
open(OUT,">",$out);
open(OUT2,">",$m);
open(IN,"<",$in);
while(<IN>){
  chomp;
  next if /^\#/;
  ($dum,$name,$chr,$str,$start,$end,$cdsStart,$cdsEnd,$dum,$exstart,$exend,$dum,$name2)=split(/\t/);
  @S = split(/\,/,$exstart);
  @E = split(/\,/,$exend);
  $sid = $name2."__".$chr;
  $t   = $name."__".$chr;
  if (!exists($count{$t})){
    $count{$t} = 1;
  }else{
    $count{$t} ++;
  }
  $tid = $name."__".$chr."__".$count{$t};
  push(@{$n{$sid}},$tid);
  if ( $utr == 0 ){
    print OUT join("\t",($chr,"UCSC","transcript",$start,$end,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
    for $i (0..$#S){
      print OUT join("\t",($chr,"UCSC","exon",$S[$i],$E[$i],0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
    }
  }else{
    if ( ( ( $utr == 5 && $str eq "+" ) || ( $utr == 3 && $str eq "-" ) )  && $start < $cdsStart ){
      print OUT join("\t",($chr,"UCSC","transcript",$start,$cdsStart - 1,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
      for $i (0..$#S){
	if ( $S[$i] < $cdsStart ){
	  if ( $E[$i] < $cdsStart ){ $e = $E[$i] }else{ $e = $cdsStart - 1};
	  print OUT join("\t",($chr,"UCSC","exon",$S[$i],$e,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
	}
      }
    }
    if ( ( ( $utr == 5 && $str eq "-" ) || ( $utr == 3 && $str eq "+" ) ) && $end > $cdsEnd ){
      print OUT join("\t",($chr,"UCSC","transcript",$cdsEnd + 1,$end,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
      for $i (0..$#S){
	if ( $E[$i] > $cdsEnd ){
	  if ( $S[$i] > $cdsEnd ){ $s = $S[$i] }else{ $s = $cdsEnd + 1};
	  print OUT join("\t",($chr,"UCSC","exon",$s,$E[$i],0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
	}
      }
    }
  }
}
close(IN);
close(OUT);

foreach ( sort keys %n ){
  print OUT2 $_."\t".join("|",@{$n{$_}})."\n";
}
close(OUT2);
