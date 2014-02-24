#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use strict;
use warnings;

my ($assemdir, $all_prot_seqs, $target_prots_list, $makeblastdb_path) = @ARGV;

# Create blast database directory
my $blast_dbs_dir = "$assemdir/blast_dbs/";
unless(-d $blast_dbs_dir or mkdir $blast_dbs_dir) {
    die "[ERROR prepareData] Could not create blast database output directory $blast_dbs_dir\n";
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
system("perl -ne 'if (/^>(\\S+)/) {\$c=\$i{\$1}}\$c?print:chomp;\$i{\$_}=1 if \@ARGV' $target_prots_list $blast_dbs_dir/all_proteins.fasta > $target_seqs");

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

=item $assemdir

Output directory.

=item $all_prot_seqs

Reference proteome file with all target protein sequences.

=item $target_prots_list

Text file with target protein IDs.

=item $makeblastdb_path

Path to makeblastdb binary.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR prepareData] Could not create blast database output directory ...

Directory creation failed. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

makeblastdb

=head1 KNOWN BUGS

None.

=head1 NOTES

makeblastdb is from the new BLAST+ package.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
