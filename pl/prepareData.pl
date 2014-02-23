# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use strict;
use warnings;

my ($assemdir, $all_prot_seqs, $target_seqs_list, $makeblastdb_path) = @ARGV;

# Create blast database directory
my $blast_dbs_dir = "$assemdir/blast_dbs/";
unless(-d $blast_dbs_dir or mkdir $blast_dbs_dir) {
    die "Could not create blast database output directory $blast_dbs_dir\n";
}

my $logfile = "$assemdir/blast_dbs/prepareData.log";

# Copy in all proteins file, using the first field as the ID
system("awk '{print \$1}' $all_prot_seqs > $blast_dbs_dir/all_proteins.fasta");

# Create the all proteins blast db
unless(-e "$blast_dbs_dir/all_proteins.pin") {
    system("$makeblastdb_path -dbtype prot -in $blast_dbs_dir/all_proteins.fasta -out $blast_dbs_dir/all_proteins > $logfile");
}

# Extract sequences in target list from all proteins file
my $target_seqs = "$blast_dbs_dir/target_proteins.fasta";
system("perl -ne 'if (/^>(\\S+)/) {\$c=\$i{\$1}}\$c?print:chomp;\$i{\$_}=1 if \@ARGV' $target_seqs_list $blast_dbs_dir/all_proteins.fasta > $target_seqs");

# Create the targets BLAST database 
unless(-e "$blast_dbs_dir/target_proteins.pin") {
    system("$makeblastdb_path -dbtype prot -in $target_seqs -out $blast_dbs_dir/target_proteins >> $logfile");
}

__END__

=head1 NAME

prepareData - Create all and target protein blast databases.

=head1 USAGE

=over

=item B<perl prepareData.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=item $all_exons_file

Nucleotide FASTA file with all target exon sequences.

=item $exonlist

Text file with target exon IDs.

=item $minoverlap

Length filtering threshold: contig-protein/exon-protein alignment range ratio must >= this value for an assembly to pass.

=item $exonerate_path

Path to exonerate binary.

=item $blastall_path

Path to blastall binary.

=back

=head1 SUBROUTINES

=over

=item getTargetRegionInProtein

Align target protein to target exon and determine alignment beginning and end indices on the protein sequence.

=item filterExoneratedContigs

Align target protein to each assembly and compare alignment range to the return values of getTargetRegionInProtein. Assemblies that fail to meet $minoverlap are discarded.

=item getBestContig

BLASTx filtered assemblies to all reference proteins; select representative sequences for target exons.

=back

=head1 DIAGNOSTICS

=over 

=item [ERROR bestcontig_distrib_dir $lib] Could not make ...

File or directory creation failed. Check that you have adequate permissions in the output directory.

=item [ERROR bestcontig_distrib_dir $lib] Could not open ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully. Check that you have adequate permissions in the output directory.

=item [WARNING bestcontig_distrib] No prot-exon exonerate alignment ...

A target protein and target exon did not align with exonerate. Check that the exon is named correctly, and is orthologous to the target protein.

=item [ERROR bestcontig_distrib_dir $lib] Invalid exonerate sequence ID line ...

Could not find alignment range in exonerate file due to malformed sequence ID line. Check for file corruption.

=item [WARNING bestcontig_distrib] No prot-contig exonerate alignment ...

Target protein did not align with exonerate to any assembled contigs. Check that previous pipeline stages completed successfully, and that there are assembled contigs in the given filename.

=back

=head1 DEPENDENCIES

exonerate

blastall

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
