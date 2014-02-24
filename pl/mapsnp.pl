#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($lib, $readdir, $assemdir, $fwd_suffix, $rev_suffix, $unpaired_suffix, $bowtie2build_path, $bowtie2_path, $samtools_path) = @ARGV;

my $mapdir = "$assemdir/$lib/${lib}_mapsnp/";
unless(-e $mapdir or mkdir $mapdir) { die "[ERROR mapsnp $lib] Unable to create $mapdir\n"; }

# bowtie2 build on representative contigs for targets
my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $bowtie2_index_name = "$mapdir/${lib}_best2refs";
system("$bowtie2build_path $reffil $bowtie2_index_name > $mapdir/${lib}_best2refs.bowtie2build.log");

my $sorted_bam = mapFiles($readdir, $lib, $bowtie2_index_name, $mapdir, $fwd_suffix, $rev_suffix, $unpaired_suffix);

# Subroutines

sub mapFiles {
    my ($readdir, $lib, $bowtie2_index_name, $mapdir, $fwd_suffix, $rev_suffix, $unpaired_suffix) = @_;

    my $file1 = $readdir . $lib . "_$fwd_suffix.fastq.gz";
    my $file2 = $readdir . $lib . "_$rev_suffix.fastq.gz";
    my $fileu = $readdir . $lib . "_$unpaired_suffix.fastq.gz";

    my $logfile = "$mapdir/$lib.mapsnp.log";
    my $sorted_bam = "$mapdir/$lib.sorted";
    system(qq(
        $bowtie2_path -5 5 -3 5 -p 1 -x $bowtie2_index_name -U $file1,$file2,$fileu 2> $logfile |
        $samtools_path view -bS - 2> $logfile | $samtools_path sort - $sorted_bam 2> $logfile
    ));
    return("$sorted_bam.bam");
}

__END__

=head1 NAME

mapsnp - Map sample reads onto best contigs with bowtie2.

=head1 USAGE

=over

=item B<perl mapsnp.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $readdir

Directory containing cleaned read files.

=item $assemdir

Output directory.

=item $fwd_suffix

Filename suffix for forward read files.

=item $rev_suffix

Filename suffix for reverse read files.

=item $unpaired_suffix

Filename suffix for unpaired read files.

=item $bowtie2build_path

Path to bowtie2-build binary.

=item $bowtie2_path

Path to bowtie2 binary.

=item $samtools_path

Path to samtools binary.

=back

=head1 SUBROUTINES

=over

=item mapFiles

Map sample reads from $readdir belonging to the sample $lib, and sort the resulting BAM file.

=back

=head1 DIAGNOSTICS

=over 

=item [ERROR mapsnp $lib] Unable to create ...

Required output location not creatable. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

bowtie2-build

bowtie2

samtools

=head1 KNOWN BUGS

None.

=head1 NOTES

Sample reads from the forward, reverse and unpaired read files are treated as unpaired due to the short reference contig lengths.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
