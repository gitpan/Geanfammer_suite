#!/usr/bin/perl
# Written by Sarah A. Teichmann on December 17th 1996.
# Use: creates a summary file from a cluster file.

$usage="\nUsage: $0 <good cluster file>";


#Create summary of cluster sizes and numbers
unless (@ARGV==1){print "$usage" && die;}
my $good_cluster_file_main=shift;

&make_clustering_summary(\$good_cluster_file_main);


#--------Subroutine: prints out summary file from cluster file.----------

#________________________________________________________________________________
# Title     : make_clustering_summary
# Usage     : &make_summ($sorted_cluster_file)
# Function  : to make a summary file of a sorted cluster file
# Example   :
# Keywords  : summary, make_cluster_summary, subclustering summary
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.4
#--------------------------------------------------------------------------------
sub make_clustering_summary{
    my ($good_cluster_file, $summary_file, @filecontent, $i, $filecontent,
        $cluster_size, @cluster_sizes, $cluster_number, $number_of_clusters,
         $summary_file, @filecontent, %hash, @keys, @temp_clu);
    $good_cluster_file=${$_[0]} || $_[0];
    $summary_file="$good_cluster_file".".summary";
    open(CLU_FILE, "$good_cluster_file");
    while(<CLU_FILE>){
          push(@temp_clu, $_);  ## copying the content to ;
          if( /^ *Cluster +size +(\d+)/i){
              $cluster_size=$1;
          }elsif (/^ *Cluster +[number]* *(\d+)/) {
              $hash{$cluster_size} ++;
          }
    }
    close(CLU_FILE);

    open(CLU_FILE, ">$good_cluster_file"); # now overwrting it.
    open (SUMM, ">$summary_file");
    print SUMM "Cluster size    No. of clusters\n";
    print CLU_FILE "Cluster size    No. of clusters\n";
    @keys=sort {$a<=>$b} keys %hash;
    for ($i=0; $i<@keys; $i++){
        print SUMM "     $keys[$i]               $hash{$keys[$i]}\n";
        print CLU_FILE "     $keys[$i]               $hash{$keys[$i]}\n";
    }
    close (SUMM);
    print CLU_FILE "\n# This file is created by $0 with make_clustering_summary sub, Details below\n\n";
    for(@temp_clu){  print CLU_FILE $_ }
    close (CLU_FILE);
    return(\$summary_file);
}

