#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;
use Bio::SeqIO;

my ($lib, $assemdir) = @ARGV;

# Open final GATK output
my $vcf = "$assemdir/$lib/${lib}_gatkSNPcalls/$lib.ReadGrouped.ReadBackedPhased.VariantAnnotated.VariantFiltered.vcf";
open(VCF, "<$vcf") or die "[ERROR vcf2ambigfasta $lib] Could not open .vcf input file $vcf\n";
my @vcf = <VCF>;
close VCF;

my $vcf2ambigfasta_dir = "$assemdir/$lib/${lib}_vcf2ambigfasta/";
unless(-e $vcf2ambigfasta_dir or mkdir $vcf2ambigfasta_dir) { die "[ERROR vcf2ambigfasta $lib] Could not make $vcf2ambigfasta_dir\n"; }

my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $ref_IN = Bio::SeqIO->new(-file => "<$reffil",
                             -format => "fasta",
                             -alphabet => "dna");

my $output_reffil = "$vcf2ambigfasta_dir/${lib}_best2refs.vcf2ambigfasta_refs.fasta";
my $ref_OUT = Bio::SeqIO->new(-file => ">$output_reffil",
                              -format => "fasta",
                              -alphabet => "dna");

# Tally types of variants
my $output_statsfil = "$vcf2ambigfasta_dir/${lib}_best2refs.vcf2ambigfasta_refs.stats";
open(STATS_OUT, ">$output_statsfil") or die "[ERROR vcf2ambigfasta $lib] Could not open output reference file $output_statsfil\n";

# Loop through reference sequences
while (my $ref = $ref_IN->next_seq) {
    my $ref_name = $ref->display_id();
    my $ref_seq = $ref->seq();

    # Whilst tallying detected variant types
    my $hets = 0;
    my $alts = 0;
    my $indel = 0;
    my $other = 0;
 
    foreach my $vcline (@vcf) {
        # a0 and a1 are the reference and alternate bases
        my ($vcf_chrom_name, $snp_pos, undef, $a0, $a1, undef, $filter_result, undef, undef, $phasing_info) = split(/\t/, $vcline);

        # Check the variant passed filtering and characterise it based on phasing information
        if ($vcf_chrom_name eq $ref_name && $filter_result eq "PASS") {
            # if indel
            if (length($a0) > 1 || length($a1) > 1) {
                $indel++;
            # if alternate homozygote, replace the reference with the alternate
            } elsif ($phasing_info =~ /^1[\|\/]1/) { 
                substr($ref_seq, $snp_pos-1, 1) = $a1;
            # if heterozygote, replace the reference with the IUPAC ambiguity code
            } elsif ($phasing_info =~ /^0[\|\/]1/) { 
                my $aa;
                if    ("$a0$a1" =~ /CT|TC/) { $aa = "Y"; }
                elsif ("$a0$a1" =~ /AG|GA/) { $aa = "R"; }
                elsif ("$a0$a1" =~ /GC|CG/) { $aa = "S"; }
                elsif ("$a0$a1" =~ /AT|TA/) { $aa = "W"; }
                elsif ("$a0$a1" =~ /GT|TG/) { $aa = "K"; }
                elsif ("$a0$a1" =~ /AC|CA/) { $aa = "M"; }
                else { die "[ERROR vcf2ambigfasta $lib] Invalid base detected in $vcline\n"}

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

__END__

=head1 NAME

vcf2ambigfasta - Write SNPs to best contigs reference.

=head1 USAGE

=over

=item B<perl vcf2ambigfasta.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR vcf2ambigfasta $lib] Could not open/make ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=item [ERROR vcf2ambigfasta $lib] Invalid base detected ...

Base combination does not have a corresponding IUPAC code. Check for .vcf file corruption.

=back

=head1 DEPENDENCIES

BioPerl

=head1 KNOWN BUGS

None.

=head1 NOTES

None.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
