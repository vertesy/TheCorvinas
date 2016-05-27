#!/usr/bin/perl -s -w

#use lib "/Users/d.grun/data/bin";
use lib "/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation";
use tools;

if (scalar @ARGV == 1)
{
    die "usage: -in=INPUT.gtf -cl=clusters_list.csv -out=OUTPUT.gtf\n" if ($ARGV[0] eq "help");
}

# What is -cl=clusters_list.csv

%stop  = ();
%strand = ();
%chr=();
open(IN,"<",$cl);
while(<IN>){
    chomp;
    ($cl_name,$e) = split(/\t/);
    @me = split(/\|/,$e);
    foreach $m (@me){
	$clust{$m} = $cl_name;
    }
}
close(IN);

open(IN,"<",$in);
while(<IN>){
    chomp;
    next if $_ =~ /^\#/;
    @F = split(/\t/);
    next if $F[2] eq "transcript";
    $info = $F[$#F];
    $info =~ s/(\"|\;)//g;
    ($dum,$dum,$dum,$tid) = split(/\s+/,$info);
    $id = $clust{$tid};
    if (!exists($stop{$id}{$F[3]})){
	$stop{$id}{$F[3]}  = $F[4];
	$strand{$id}       = $F[6];
	$chr{$id}          = $F[0];
    }else{
	$stop{$id}{$F[3]} = $F[4] if $F[4] > $stop{$id}{$F[3]};
    }
}
close(IN);
open(OUT,">",$out);
foreach $id (sort keys %stop){
    @start    = sort {$a <=> $b} keys %{$stop{$id}};
    $tmp_stop = -1;
    @exon     = ();

    foreach $i (0..$#start){
	$s = $start[$i];
	if ($i == 0){
	    $tmp_min = $s;
	    $tmp_max = $stop{$id}{$s};
	    $tmp_start = $s;
	    $tmp_stop = $stop{$id}{$s};
	}else{
	    $tmp_min = min($tmp_min, $s);
	    $tmp_max = max($tmp_max,$stop{$id}{$s});
	    if ($s > $tmp_stop ){
		push(@exon,join("\t",($chr{$id},"merge","exon",$tmp_start,$tmp_stop,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\"; exon \"".$i."\";"))."\n");
		$tmp_start = $s;
	    }
	    $tmp_stop  = max( $tmp_stop, $stop{$id}{$s});
	}
    }
    push(@exon,join("\t",($chr{$id},"merge","exon",$tmp_start,$tmp_stop,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\"; exon \"".($#start + 1)."\";"))."\n");
    print OUT join("\t",($chr{$id},"merge","transcript",$tmp_min,$tmp_max,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\""))."\n";
    foreach $e (@exon){
	print OUT $e;
    }
}
close(OUT);
