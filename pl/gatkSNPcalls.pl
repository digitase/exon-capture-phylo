use warnings;
use strict;

# my $readdir  = "/home2/jgb/reads/";
# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $mapdir   = "/home2/jgb/assemble/crypto/map/";
# my $lib = $ARGV[0];

my ($lib, $readdir, $assemdir, $picard_dir, $gatk_dir, $java_heap_size, $samtools_path) = @ARGV;

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
    my $logfile = "$gatkSNPcalls_dir/$lib.picard.log";

    # Add read groups to BAM
    # Set max heap size
    my $AddOrRepl = "java -Xmx${java_heap_size}g -jar $picard_dir/AddOrReplaceReadGroups.jar";
    system("$AddOrRepl INPUT=$bam OUTPUT=$ibamrg RGID=$lib RGLB=$lib RGPU=$lane RGPL=illumina RGSM=$samp 2> $logfile");   

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

    # Variant calls
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

    # DepthOfCoverage
    system("$gatk -R $ref -T DepthOfCoverage -I $ibamrg -o $ibamrg_basepath.DepthOfCoverageTable >> $logfile"); 

    return($filtered_vcf);
}

__END__;



