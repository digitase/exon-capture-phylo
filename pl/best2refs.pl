#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($lib, $assemdir, $exonlist) = @ARGV;

my $best2refs_dir = "$assemdir/$lib/${lib}_best2refs/";
unless(-e $best2refs_dir or mkdir $best2refs_dir) { die "[ERROR best2refs $lib] Could not make $best2refs_dir\n"; }

# Open target exons list
open EXONS, "<$exonlist" or die "[ERROR best2refs $lib] Could not open the exons list file $exonlist\n";
my @exons = <EXONS>;
close(EXONS);
chomp(@exons);

my $gathercontigs_dir = "$assemdir/gathercontigs_by_target_exon/";

# Open output file
my $reffil = "$best2refs_dir/${lib}_best2refs.fasta";
open REFS, ">$reffil" or die "[ERROR best2refs $lib] Could not open the output best2refs fasta file $reffil\n";

foreach my $exon_name (@exons) {

    if ($exon_name =~ /^(\S+?)_\S+/) { 

        # Non-whitespace portion of exon ID before underscore should be the orthologous protein ID
        my $prot = $1; 

        # Get representative contig for target exon
        my $bestcontig_distrib_dir = "$assemdir/$lib/$prot/${prot}_bestcontig_distrib/";
        my $bestcontigfil = "$bestcontig_distrib_dir/${exon_name}_velvet_contigs.cap3ed.exonerated.filtered.best_contig.fasta";      

        # Read contig into array
        open CTG, "<$bestcontigfil" or die "[ERROR best2refs $lib] Could not open the best contig file at $bestcontigfil\n";
        my @contig = <CTG>;
        chomp(@contig);
        close(CTG); 

        # Rename the sequence
        if (scalar(@contig) == 2) {
            my $seq = $contig[1];
            print REFS ">$exon_name\n" . "$seq\n"; 
        } elsif (scalar(@contig) == 0) {
            # This is fine. There may be no representative sequence due to failed assembly or filtering.
            warn "[WARNING best2refs $lib] No best contig for $exon_name\n";
        } else {
            die "[ERROR best2refs $lib] Incorrect number of contigs/contig format in $bestcontigfil\n";
        }

    } else {
        die "[ERROR best2refs $lib] Invalid exon name $exon_name found in $exonlist\n";
    }
    
}

__END__

=head1 NAME

best2refs - Gather representative sequences for each target exon into one file, for use as a reference during SNP calling.

=head1 USAGE

=over

=item B<perl best2refs.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=item $exonlist

Text file with target exon IDs.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR best2refs $lib] Could not open/make ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=item [WARNING best2refs $lib] No best contig ...

No representative sequence located. May have failed to assemble, or was filtered out at the bestcontig_distrib.pl stage.

=item [ERROR best2refs $lib] Incorrect number of contigs/contig format ...

FASTA record of best contig is not in the correct format. Check for file corruption.

=item [ERROR best2refs $lib] Invalid exon name ...

Exon naming format does not conform to specifications.

=back

=head1 DEPENDENCIES

None.

=head1 KNOWN BUGS

None.

=head1 NOTES

The regex used to parse target exon names restricts the sequence ID format:
    if ($exon_name =~ /^(\S+?)_\S+/) 

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
