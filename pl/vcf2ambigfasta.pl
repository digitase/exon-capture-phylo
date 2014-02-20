#!/usr/bin/perl -w

use warnings;
use strict;
use Bio::SeqIO;

my ($lib, $assemdir) = @ARGV;

my $vcf = "$assemdir/$lib/${lib}_gatkSNPcalls/$lib.ReadGrouped.ReadBackedPhased.VariantAnnotated.VariantFiltered.vcf";
open(VCF, "<$vcf") or die "Could not open .vcf input file $vcf\n";
my @vcf = <VCF>;
close VCF;

my $vcf2ambigfasta_dir = "$assemdir/$lib/${lib}_vcf2ambigfasta/";
unless(-e $vcf2ambigfasta_dir or mkdir $vcf2ambigfasta_dir) { die "Unable to create $vcf2ambigfasta_dir\n"; }

my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $ref_IN = Bio::SeqIO->new(-file => "<$reffil",
                             -format => "fasta",
                             -alphabet => "dna");

my $output_reffil = "$vcf2ambigfasta_dir/${lib}_best2refs.vcf2ambigfasta_refs.fasta";
my $ref_OUT = Bio::SeqIO->new(-file => ">$output_reffil",
                              -format => "fasta",
                              -alphabet => "dna");

my $output_statsfil = "$vcf2ambigfasta_dir/${lib}_best2refs.vcf2ambigfasta_refs.stats";
open(STATS_OUT, ">$output_statsfil") or die "Could not open output reference file $output_statsfil\n";

# loop through reference sequences
while (my $ref = $ref_IN->next_seq) {
    my $ref_name = $ref->display_id();
    my $ref_seq = $ref->seq();

    my $hets = 0;
    my $alts = 0;
    my $indel = 0;
    my $other = 0;
 
# 0                           1   2   3   4   5       6       7       
# CHROM                       POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  indexing4
# ENSACAP00000021572_exon1    732 .   C   T   2344.77 PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-1.664;ClippingRankSum=-0.732;DP=153;FS=16.321;HaplotypeScore=6.6471;MLEAC=1;MLEAF=0.500;MQ=41.67;MQ0=0;MQRankSum=1.320;QD=15.33;ReadPosRankSum=1.222   GT:AD:GQ:PL:PQ  1|0:96,89:99:2373,0,2303:2115.71
    foreach my $vcline (@vcf) {
        my ($vcf_chrom_name, $snp_pos, undef, $a0, $a1, undef, $filter_result, undef, undef, $phasing_info) = split(/\t/, $vcline);

        # check the snp fits criteria
        if ($vcf_chrom_name eq $ref_name && $filter_result eq "PASS") {
            # indel
            if (length($a0) > 1 || length($a1) > 1) {
                $indel++;
            # if alternate homozygote         
            } elsif ($phasing_info =~ /^1[\|\/]1/) { 
                substr($ref_seq, $snp_pos-1, 1) = $a1;
            # if heterozygote             
            } elsif ($phasing_info =~ /^0[\|\/]1/) { 
                my $aa;

                if    ("$a0$a1" =~ /CT|TC/) { $aa = "Y"; }
                elsif ("$a0$a1" =~ /AG|GA/) { $aa = "R"; }
                elsif ("$a0$a1" =~ /GC|CG/) { $aa = "S"; }
                elsif ("$a0$a1" =~ /AT|TA/) { $aa = "W"; }
                elsif ("$a0$a1" =~ /GT|TG/) { $aa = "K"; }
                elsif ("$a0$a1" =~ /AC|CA/) { $aa = "M"; }
                else { die "Invalid base detected in $vcline\n"}

                substr($ref_seq, $snp_pos-1, 1) = $aa;
                $hets++;
            } else {
                $other++;
            }
        }
    }

    $ref->seq("$ref_seq");
    my $seqlen = $ref->length();
    $ref_OUT->width($seqlen);
    $ref_OUT->write_seq($ref);

    print STATS_OUT "$ref_name\tHET:$hets\tALT_HOMO:$alts\tINDEL:$indel\tSEQ_LEN:$seqlen\n";
}

close STATS_OUT;
