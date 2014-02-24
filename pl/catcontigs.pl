#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($lib, $assemdir, $target_seqs_list, $cap3_path, $exonerate_path, @k_values) = @ARGV;

my $assemlib = "$assemdir/$lib/";

# Read list of target proteins
open PROTS, "<$target_seqs_list" or die "[ERROR catcontigs $lib] Could not open target proteins list $target_seqs_list\n";
my @prots = <PROTS>;
chomp(@prots); 
close(PROTS);

foreach my $prot (@prots) {

    my $catcontigs_dir = "$assemlib/$prot/${prot}_catcontigs/";
    unless(-e $catcontigs_dir or mkdir $catcontigs_dir) { die "[ERROR catcontigs $lib] Could not make $catcontigs_dir\n"; }

    # Copy in velvet contigs from all k values
    my $all_assembled_contigs = "$assemlib/$prot/${prot}_call_velvet_assemblies/${prot}_velvet_contigs.fasta";
    system("cp $all_assembled_contigs $catcontigs_dir");
    $all_assembled_contigs = "$catcontigs_dir/${prot}_velvet_contigs.fasta";

    # Require 99% identity in a 20 bp overlap to merge with CAP3
    system("$cap3_path $all_assembled_contigs -o 20 -p 99 > $catcontigs_dir/${prot}_cap3.log");

    # Collate consensus and singleton contigs from CAP3
    my $cap3ed_contigs = "$catcontigs_dir/${prot}_velvet_contigs.cap3ed.fasta";
    system("cat $all_assembled_contigs.cap.contigs $all_assembled_contigs.cap.singlets > $cap3ed_contigs");

    # Align target protein to assembled contigs with exonerate
    my $exonerate_query = "$catcontigs_dir/$prot.fasta";
    system("perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($prot) } print if \$c' $assemdir/blast_dbs/target_proteins.fasta > $exonerate_query");

    # Define custom exonerate output format that includes alignment range
    my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
    my $exonerated_contigs = "$catcontigs_dir/${prot}_velvet_contigs.cap3ed.exonerated.fasta";
    system("$exonerate_path --model protein2genome --query $exonerate_query --target $cap3ed_contigs --showvulgar no --showalignment no --ryo $ryo > $exonerated_contigs");

    # Uncomment to reduce IO demands
    # sleep(2);
}

__END__

=head1 NAME

catcontigs - Merge assembled contigs with CAP3. Align assembled contigs to target protein with exonerate and output alignment range. 

=head1 USAGE

=over

=item B<perl catcontigs.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=item $target_seqs_list

Text file with target protein IDs.

=item $cap3_path

Path to cap3 binary.

=item $exonerate_path

Path to exonerate binary.

=item $k_values

Array of Velvet k-values that were used. Should be given as a Bash array e.g.
    VELVET_K_VALUES=(31 41 51 61)
    perl catcontigs.pl [Other ARGS] ${VELVET_K_VALUES[@]}

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR catcontigs $lib] Could not make ...

File or directory creation failed. Check that you have adequate permissions in the output directory.

=item [ERROR catcontigs $lib] Could not open target proteins list ...

Check that TARGET_PROTEIN_SEQS_LIST is set correctly in the config file.

=back

=head1 DEPENDENCIES

cap3

exonerate

=head1 KNOWN BUGS

None.

=head1 NOTES

Exonerate alignment range is reported as the start and end of the alignment region in the query sequence (the target protein). For more information on exonerate "ryo" formats, see https://www.ebi.ac.uk/~guy/exonerate/advanced.html and http://csurs7.csr.uky.edu/cgi-bin/man/man2html?1+exonerate

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
