#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($assemdir, $samples_list, $exonlist) = @ARGV;

my $gatherambigcontigs_dir = "$assemdir/gatherAmbigcontigs/";
unless(-e $gatherambigcontigs_dir or mkdir $gatherambigcontigs_dir) { die "[ERROR gatherambigcontigs] Could not make $gatherambigcontigs_dir\n"; }    


open LIBS, "<$samples_list" or die "[ERROR gatherAmbigcontigs] Could not open the sample names file $samples_list\n";
my @libs  = <LIBS>;
close(LIBS);
chomp(@libs);

open EXONS, "<$exonlist" or die "[ERROR gatherAmbigcontigs] Could not open the target exon IDs file $exonlist\n";
my @exons = <EXONS>;
close(EXONS);
chomp(@exons);

foreach my $exon_name (@exons) {
    
    # Output file for each target exon
    my $exon_best_ambig_contigs = "$gatherambigcontigs_dir/${exon_name}_best_ambig_contigs.fasta";

    foreach my $lib (@libs) {

        # For each sample, extract best contig sequence and change sequence ID to sample name
        my $vcf2ambigfasta_refs_fil = "$assemdir/$lib/${lib}_vcf2ambigfasta/${lib}_best2refs.vcf2ambigfasta_refs.fasta";
        system(qq(
            perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($exon_name) } print if \$c' $vcf2ambigfasta_refs_fil |
            sed 's/^>.*/>$lib/' >> $exon_best_ambig_contigs
        ));
    }
}

__END__

=head1 NAME

gatherAmbigcontigs - For each target, from each sample, gather the best contig sequence that was amended with IUPAC codes in vcf2ambigfasta.pl

=head1 USAGE

=over

=item B<perl gathercontigs.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $assemdir

Output directory.

=item $samples_list

List of sample names.

=item $exonlist

Text file with target exon IDs.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR gathercontigs] Could not open/make ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

None.

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
