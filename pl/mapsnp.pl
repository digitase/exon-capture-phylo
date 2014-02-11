use warnings;
use strict;

# my $readdir  = "/home2/jgb/reads/";
# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $lib = $ARGV[0];

my ($lib, $readdir, $assemdir, $bowtie2build_path, $bowtie2_path, $samtools_path) = @ARGV;

my $mapdir = "$assemdir/$lib/${lib}_mapsnp/";
unless(-e $mapdir or mkdir $mapdir) { die "Unable to create $mapdir\n"; }

my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
my $bowtie2_index_name = "$mapdir/${lib}_best2refs";
system("$bowtie2build_path $reffil $bowtie2_index_name > $mapdir/${lib}_best2refs.bowtie2build.log");

my $sorted_bam = mapFiles($readdir, $lib, $bowtie2_index_name, $mapdir);

# TODO unhardcode file formats
sub mapFiles {
    my ($readdir, $lib, $bowtie2_index_name, $mapdir) = @_;

    my $file1 = $readdir . $lib . '_1_final.fastq.gz';
    my $file2 = $file1; $file2 =~ s/_1_/_2_/;
    my $fileu = $file1; $fileu =~ s/_1_/_u_/;

    my $logfile = "$mapdir/$lib.mapsnp.log";
    my $sorted_bam = "$mapdir/$lib.sorted";
    # Treat reads as unpaired due to small references
    system(qq(
        $bowtie2_path -5 5 -3 5 -p 1 -x $bowtie2_index_name -U $file1,$file2,$fileu 2> $logfile |
        $samtools_path view -bS - 2> $logfile | $samtools_path sort - $sorted_bam 2> $logfile
    ));
    return("$sorted_bam.bam");
}

