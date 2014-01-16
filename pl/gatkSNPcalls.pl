use warnings;
use strict;

# my $readdir  = "/home2/jgb/reads/";
# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $mapdir   = "/home2/jgb/assemble/crypto/map/";
# my $lib = $ARGV[0];

my ($lib, $readdir, $assemdir, $picard_dir, $gatk_dir) = @ARGV;

my $gatkSNPcalls_dir = "$assemdir/$lib/${lib}_gatkSNPcalls/";
unless(-e $gatkSNPcalls_dir or mkdir $gatkSNPcalls_dir) { die "Unable to create $gatkSNPcalls_dir\n"; }

my $mapdir = "$assemdir/$lib/${lib}_mapsnp/";
my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $sorted_bam = "$mapdir/$lib.sorted.bam";

my $bamrg  = prepareBAMandRef($lib, $reffil, $sorted_bam, $gatkSNPcalls_dir, $picard_dir);
my $vcffilt = callGATK($bamrg, $reffil, $gatk_dir, $gatkSNPcalls_dir, $lib);

sub prepareBAMandRef {
    my ($lib, $ref, $bam, $gatkSNPcalls_dir, $picard_dir) = @_;

    # Sample names are of this format SP04_indexing12
    my ($lane, $samp) = split(/_/, $lib);
    my $ibamrg = "$gatkSNPcalls_dir/$lib.ReadGrouped.bam";

    # Add read groups to BAM
    # Set max heap size 8G
    my $AddOrRepl = "java -Xmx8g -jar $picard_dir/AddOrReplaceReadGroups.jar";
    system("$AddOrRepl INPUT=$bam OUTPUT=$ibamrg RGID=$lib RGLB=$lib RGPU=$lane RGPL=illumina RGSM=$samp");   

    # Index BAM file
    system("samtools index $ibamrg");

    # Creating the fasta sequence dictionary file
    (my $dict = $ref) =~ s/\.fasta$/\.dict/;
    unless (-e $dict) { system("java -Xmx8g -jar $picard_dir/CreateSequenceDictionary.jar R=$ref O=$dict"); }

    # Index reference FASTA
    system("samtools faidx $ref");

    return($ibamrg);
}

sub callGATK {
    my ($ibamrg, $ref, $gatk_dir, $gatkSNPcalls_dir, $lib) = @_;

    my $gatk  = "java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar";

    # Variant calls
    my $ibamrg_basepath  = $ibamrg;
    $ibamrg_basepath =~ s/\.bam$//;

    my $vcf = "$ibamrg_basepath.vcf";
    system("$gatk -R $ref -T HaplotypeCaller  -I $ibamrg -o $vcf");

    my $phased_vcf = "$ibamrg_basepath.ReadBackedPhased.vcf";
    system("$gatk -R $ref -T ReadBackedPhasing -I $ibamrg --variant $vcf --min_base_quality_score 21 -o $phased_vcf");

    my $annotated_vcf = "$ibamrg_basepath.ReadBackedPhased.VariantAnnotated.vcf";
    system("$gatk -R $ref -T VariantAnnotator  -I $ibamrg -A DepthPerAlleleBySample -A HaplotypeScore --variant $phased_vcf -o $annotated_vcf");

    my $filtered_vcf = "$ibamrg_basepath.ReadBackedPhased.VariantAnnotated.VariantFiltered.vcf";
    system("$gatk -R $ref -T VariantFiltration  -o $filtered_vcf --variant $annotated_vcf --filterName depth --filterExpression \"DP \< 16\"");

    # DepthOfCoverage
    system("$gatk -R $ref -T DepthOfCoverage -I $ibamrg -o $ibamrg_basepath.DepthOfCoverageTable"); 

    return($filtered_vcf);
}

__END__;



