use warnings;
use strict;

my $readdir  = "/home2/jgb/reads/";
my $assemdir = "/home2/jgb/assemble/crypto/";
my $mapdir   = "/home2/jgb/assemble/crypto/map/";
my $lib = $ARGV[0];
my $reffil = $assemdir . $lib . "/" . $lib . "_refs.fasta";
 
my $refs    = bw($reffil);
my $bam     = mapFiles($readdir,$lib,$refs,$mapdir);
callVariants($lib,$mapdir,$bam,$refs);




sub callVariants {
    my ($lib,$mapdir,$bam,$refs) = @_;
    
    my $outVCF = $mapdir . $lib . '.vcf';

    my $call1 = system("samtools faidx $refs");
    #my $call2 = system("samtools mpileup -uf $refs $bam | bcftools view -bvcg - > $lib.var.raw.bcf");
    #my $call3 = system("bcftools view $lib.var.raw.bcf | vcfutils.pl varFilter -w 0 > $outVCF");
}

sub mapFiles {
    my ($readdir,$lib,$refs,$mapdir) = @_;
    my $file1 = $readdir . $lib . '_1_final.fastq.gz';
    my $file2 = $file1; $file2 =~ s/_1_/_2_/;
    my $fileu = $file1; $fileu =~ s/_1_/_u_/;

    my $sam      = $mapdir . $lib . ".sam";
    my $bam      = $mapdir . $lib . ".bam";
    my $finalout = $mapdir . $lib . ".sorted";

    my $call1 = system("bowtie2 -5 5 -3 5 -p 1 -x $refs -U $file1,$file2,$fileu -S $sam");
    my $call2 = system("samtools view -bS $sam > $bam");
    my $call3 = system("samtools sort $bam $finalout");
    $finalout .= ".bam";
    return($finalout);
}

sub bw {
    my $file = $_[0];
    unless(-f $file . ".1.bt2") {
    my $call = system("/usr/local/bin/bowtie2-build $file $file");
    }
    return($file);
}
