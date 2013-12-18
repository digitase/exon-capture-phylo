#!/usr/bin/perl -w

use warnings;
use strict;
use Bio::SeqIO;


my $assemdir = "/home2/jgb/assemble/target/";
my $mapdir   = "/home2/jgb/assemble/target/map/";

my $lib = $ARGV[0];

my $reffil = $assemdir . $lib . "/" . $lib . "_refs.fasta";
my $seqfil = $assemdir . $lib . "/" . $lib . "_seqs.fasta";
my $vcf    = $assemdir . $lib . "/" . $lib . ".sortedgatk.vcf.filt";
my $vcfsnp = $vcf . ".snp";

my $statfil= $seqfil . ".stat";

my %loci;

open(SEQ, ">$seqfil");
open(STAT, ">$statfil");

#my $gatk  = "java -Xmx8g -jar /home/jgb/software/GenomeAnalysisTK-nightly-2013-04-11-gb82c674/GenomeAnalysisTK.jar";
#my $sv = "SelectVariants";
#system("$gatk -R $ref -T $sv -o $vcfsnp --variant $vcf -selectType SNP -restrictAllelesTo BIALLELIC");
#open(VCF, "<$vcfsnp");

open(VCF, "<$vcf");
my @vcf = <VCF>;
 
  # handle IO with bioperl
  my $seqio = Bio::SeqIO-> new(
                               -file     => $ref,
                               -format => 'FASTA',
                               );

  # loop through reference sequences
  while( my $seq = $seqio->next_seq ) {
  
    my $name = $seq->id;
    my $nameseq = $seq->seq();

    my $hets = 0;
    my $alts = 0;
    my $other = 0;
 
    # get vcs for contig
    my @refvcs = ();

    foreach my $vcline (@vcf)   
    {
       my @tmp = split(/\t/,$vcline);
       if ($tmp[6] eq "PASS" )  # check the snp fits criteria
       {
 
          my $a0; my $a1; my $aa;      
          if ($name eq $tmp[0])
          {              
            $a0 = $tmp[3];
            $a1 = $tmp[4];
            
            if ($tmp[9] ~= /^0[\|\/]1/) # if heterozygote             
            { 
              if (( $a0 eq "C" && $a1 eq "T" ) || ( $a1 eq "C" && $a0 eq "T" )) { $aa = "Y"; }
              if (( $a0 eq "A" && $a1 eq "G" ) || ( $a1 eq "A" && $a0 eq "G" )) { $aa = "R"; }
              if (( $a0 eq "A" && $a1 eq "T" ) || ( $a1 eq "A" && $a0 eq "T" )) { $aa = "W"; }
              if (( $a0 eq "G" && $a1 eq "C" ) || ( $a1 eq "G" && $a0 eq "C" )) { $aa = "S"; }
              if (( $a0 eq "T" && $a1 eq "G" ) || ( $a1 eq "T" && $a0 eq "G" )) { $aa = "K"; }
              if (( $a0 eq "C" && $a1 eq "A" ) || ( $a1 eq "C" && $a0 eq "A" )) { $aa = "M"; }
              $hets++
            }

            if ($tmp[9] ~= /^1[\|\/]1/) # if alternate homozygote         
            { 
              if ( $a1 eq "A" || $a1 eq "C" || || $a1 eq "G" || || $a1 eq "T" ) {$aa = $a1;}
              $alts++
            }

            if ( length($a0) > 1 || length($a1) > 1 )
            { $other++ }

            substr($nameseq, $tmp[1]-1, 1) = $aa;

         } 
       }
    }

    print SEQ "> " . $name . "\n" . $nameseq . "\n";
    my $nameseqlen = length($nameseq);
    print STAT $name . "\t" . $hets . "\t" . $alts . "\t" . $other . "\t" . $nameseqlen . "\n";
  }
  close SEQ;
  close STAT;
}
close VCF;
