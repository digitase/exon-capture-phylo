# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($lib, $assemdir, $target_seqs_list, $velveth_path, $velvetg_path, $k) = @ARGV;

my $time = localtime();
print "callVelvetAssemblies on $lib with k = $k at $time\n";

my $assemlib = "$assemdir/$lib/";

# Read list of target proteins
open TARGET_SEQS_LIST, "<$target_seqs_list" or die "[ERROR callVelvetAssemblies $lib $k] Could not open target protein IDs list $target_seqs_list";
my @protnames = <TARGET_SEQS_LIST>;
chomp(@protnames);
close(TARGET_SEQS_LIST);

foreach my $protname (@protnames) {
    # Create directory for specific k-value
    my $veldir = "$assemlib/$protname/${protname}_call_velvet_assemblies/";
    unless(-e $veldir or mkdir $veldir) { die "[ERROR callVelvetAssemblies $lib $k] Could not make $veldir\n"; }

    my $kvalue_dir = "$veldir/${protname}_k$k/";
    unless(-e $kvalue_dir or mkdir $kvalue_dir) { die "[ERROR callVelvetAssemblies $lib $k] Could not make $kvalue_dir\n"; }

    # Collate hit reads originating from all sample read files
    my $allhitreads  = "$kvalue_dir/${protname}_all_hitreads.fasta";
    system("cat $assemlib/$protname/${protname}_assemble_by_prot/* > $allhitreads");

    # Assemble with parameters optimised for small assemblies
    system("$velveth_path $kvalue_dir $k -short -fasta $allhitreads > $kvalue_dir/${protname}_velveth.log");
    system("$velvetg_path $kvalue_dir -very_clean yes -max_branch_length 320 -max_gap_count 6 -cov_cutoff 5 > $kvalue_dir/${protname}_velvetg.log");

    # Append assemblies to file that collates assemblies for all k-values
    my $all_assembled_contigs = "$veldir/${protname}_velvet_contigs.fasta";
    system("cat $kvalue_dir/contigs.fa >> $all_assembled_contigs");

    # Uncomment to reduce IO demands
    # sleep(2); 
}

__END__

=head1 NAME

callVelvetAssemblies - Assemble sample reads aligned to each target protein at a specific k-value with Velvet. 

=head1 USAGE

=over

=item B<perl callVelvetAssemblies.pl [ARGS]>

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

=item $velveth_path

Path to velveth binary.

=item $velvetg_path

Path to velvetg binary.

=item $k

Velvet k-value to use.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR callVelvetAssemblies $lib $k] ... 

File or directory creation failed. Check that you have adequate permissions in the output directory.

=item [ERROR callVelvetAssemblies $lib $k] Could not open target protein IDs list...

Check that TARGET_PROTEIN_SEQS_LIST is set correctly in the config file.

=back

=head1 DEPENDENCIES

velveth

velvetg

=head1 KNOWN BUGS

None.

=head1 NOTES

This script assembles at one k-value at a time; a for loop or xargs can be used to iterate over the range of k-values.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
