# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($lib, $readdir, $assemdir, $picard_dir, $gatk_dir, $java_heap_size, $samtools_path) = @ARGV;

my $gatkSNPcalls_dir = "$assemdir/$lib/${lib}_gatkSNPcalls/";
unless(-e $gatkSNPcalls_dir or mkdir $gatkSNPcalls_dir) { die "[ERROR gatkSNPcalls $lib] Unable to create $gatkSNPcalls_dir\n"; }

# Locations of previous pipeline outputs
my $mapdir = "$assemdir/$lib/${lib}_mapsnp/";
my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $sorted_bam = "$mapdir/$lib.sorted.bam";

# Prepare inputs for GATK
my $bamrg  = prepareBAMandRef($lib, $reffil, $sorted_bam, $gatkSNPcalls_dir, $picard_dir);

# Use GATK to call variants
my $vcffilt = callGATK($bamrg, $reffil, $gatk_dir, $gatkSNPcalls_dir, $lib);

# Subroutines

sub prepareBAMandRef {
    my ($lib, $ref, $bam, $gatkSNPcalls_dir, $picard_dir) = @_;

    my $ibamrg = "$gatkSNPcalls_dir/$lib.ReadGrouped.bam";
    my $logfile = "$gatkSNPcalls_dir/$lib.picard.log";

    # Add a read group to BAM
    my $AddOrRepl = "java -Xmx${java_heap_size}g -jar $picard_dir/AddOrReplaceReadGroups.jar";
    system("$AddOrRepl INPUT=$bam OUTPUT=$ibamrg RGID=$lib RGLB=$lib RGPU=$lib RGPL=illumina RGSM=$lib 2> $logfile");   

    # Index BAM file
    system("$samtools_path index $ibamrg");

    # Creating the fasta sequence dictionary file
    (my $dict = $ref) =~ s/\.fasta$/\.dict/;
    unless (-e $dict) { system("java -Xmx${java_heap_size}g -jar $picard_dir/CreateSequenceDictionary.jar R=$ref O=$dict 2> $logfile"); }

    # Index reference FASTA
    system("$samtools_path faidx $ref");

    return($ibamrg);
}

sub callGATK {
    my ($ibamrg, $ref, $gatk_dir, $gatkSNPcalls_dir, $lib) = @_;

    my $gatk  = "java -Xmx${java_heap_size}g -jar $gatk_dir/GenomeAnalysisTK.jar";

    # Variant calling pipeline:
    my $ibamrg_basepath  = $ibamrg;
    $ibamrg_basepath =~ s/\.bam$//;

    my $logfile = "$gatkSNPcalls_dir/$lib.gatk.log";

    my $vcf = "$ibamrg_basepath.vcf";
    system("$gatk -R $ref -T HaplotypeCaller  -I $ibamrg -o $vcf >> $logfile");

    my $phased_vcf = "$ibamrg_basepath.ReadBackedPhased.vcf";
    system("$gatk -R $ref -T ReadBackedPhasing -I $ibamrg --variant $vcf --min_base_quality_score 21 -o $phased_vcf >> $logfile");

    my $annotated_vcf = "$ibamrg_basepath.ReadBackedPhased.VariantAnnotated.vcf";
    system("$gatk -R $ref -T VariantAnnotator  -I $ibamrg -A DepthPerAlleleBySample -A HaplotypeScore --variant $phased_vcf -o $annotated_vcf >> $logfile");

    my $filtered_vcf = "$ibamrg_basepath.ReadBackedPhased.VariantAnnotated.VariantFiltered.vcf";
    system("$gatk -R $ref -T VariantFiltration  -o $filtered_vcf --variant $annotated_vcf --filterName depth --filterExpression \"DP \< 16\" >> $logfile");

    # DepthOfCoverage analysis
    system("$gatk -R $ref -T DepthOfCoverage -I $ibamrg -o $ibamrg_basepath.DepthOfCoverageTable >> $logfile"); 

    return($filtered_vcf);
}

__END__

=head1 NAME

gatkSNPcalls - Call, annotate and filter variants with GATK.

=head1 USAGE

=over

=item B<perl gatkSNPcalls.pl [ARGS]>

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

=item $picard_dir 

Path to Picard tools .jar directory.

=item $gatk_dir 

Path to GATK .jar directory.

=item $java_heap_size

Heap size (gigabytes) to allocate to the Java VM for Picard and GATK.

=item $samtools_path

Path to samtools binary.

=back

=head1 SUBROUTINES

=over

=item prepareBAMandRef

Use Picard and samtools to read group and index the BAM file, and index and create a sequence dictionary from the best contigs reference file.

=item callGATK

Use GATK to call, phase, annotate and filter variants; and analyse mapping coverage in the BAM file.

=back

=head1 DIAGNOSTICS

=over 

=item [ERROR gatkSNPcalls $lib] Unable to create ...

Required output location not creatable. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

Picard

samtools

GATK

=head1 KNOWN BUGS

None.

=head1 NOTES

None.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

Refer to http://www.broadinstitute.org/gatk/gatkdocs/ for documentation on individual GATK operations.

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
