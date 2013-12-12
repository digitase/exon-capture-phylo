use warnings;
use strict;

my $readdir  = "/home2/jgb/reads/";
my $assemdir = "/home2/jgb/assemble/crypto/";
my $mapdir   = "/home2/jgb/assemble/crypto/map/";

my $lib = $ARGV[0];

my $reffil = $assemdir . $lib . "/" . $lib . "_refs.fasta";

my $bam    = $mapdir . $lib . ".sorted.bam";
my $bamrg  = rg($lib, $bam);
my $ibamrg = indbam($bamrg);

my $vcffilt = callGATK($ibamrg, $reffil);



# java -jar /home/jgb/software/picard/picard-tools-1.88/CreateSequenceDictionary.jar R=SP02B_indexing26_finalref.fa O=SP02B_indexing26_finalref.dict
# samtools faidx SP02B_indexing26_finalref.fa

sub rg {
    my ($lib, $bam)   = @_;
    my $bamrg         = $bam . ".rg";
    my ($lane, $samp) = split(/_/,$lib);

    my $AddOrRepl = "java -Xmx8g -jar /home/jgb/software/picard-tools-1.88/AddOrReplaceReadGroups.jar";
    system("$AddOrRepl INPUT=$bam OUTPUT=$bamrg RGID=$lib RGLB=$lib RGPU=$lane RGPL=illumina RGSM=$samp");   
    return($bamrg);
}


sub indbam {
    my $bamrg   = $_[0];
    my $ibamrg = $bamrg;
    $ibamrg =~ s/\.bam\.rg//;
    $ibamrg = $ibamrg . "gatk.bam";
    system("mv $bamrg $ibamrg");
    system("samtools index $ibamrg");
    return($ibamrg);
}



sub callGATK {
    my ($ibamrg, $ref) = @_;
    my $seqdict = "java -Xmx8g -jar /home/jgb/software/picard-tools-1.88/CreateSequenceDictionary.jar";
    my $gatk  = "java -Xmx8g -jar /home/jgb/software/GenomeAnalysisTK-nightly-2013-04-11-gb82c674/GenomeAnalysisTK.jar";
    my $dict  = $ref;
       $dict  =~ s/\.fasta/\.dict/; 
    my $ug    = "UnifiedGenotyper";
    my $hc    = "HaplotypeCaller";
    my $rbp   = "ReadBackedPhasing";
    my $va    = "VariantAnnotator";
    my $vf    = "VariantFiltration";
    my $doc   = "DepthOfCoverage";

    my $vcf  = $ibamrg;
    $vcf =~ s/\.bam/\.vcf/;

    unless (-e $dict) { system("$seqdict R=$ref O=$dict"); }
    system("$gatk -R $ref -T $hc  -I $ibamrg -o $vcf");
    system("$gatk -R $ref -T $rbp -I $ibamrg --variant $vcf --min_base_quality_score 21 -o $vcf.rbp");
    system("$gatk -R $ref -T $va  -I $ibamrg -A DepthPerAlleleBySample -A HaplotypeScore --variant $vcf.rbp -o $vcf.annot");
    system("$gatk -R $ref -T $vf  -o $vcf.filt --variant $vcf.annot --filterName depth --filterExpression \"DP \< 16\"");
    system("$gatk -R $ref -T $doc -I $ibamrg -o $ibamrg.cov"); 

    my $vcffilt = "$vcf.filt";

    return($vcffilt);
}



__END__;



