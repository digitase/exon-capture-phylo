use warnings;
use strict;

# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $libfil   = $assemdir . "libraries.txt";
# my $exondir  = $assemdir . "refs/targets/";
# my $exonlist = $assemdir . "refs/targets/targetexons.txt.all";

my ($assemdir, $libfil, $exondir, $exonlist) = @ARGV;

open LIBS, "<$libfil" or die "could not open the lib file";
open EXONS, "<$exonlist" or die "could not open the lib file";

my @libs  = <LIBS>;
my @exons = <EXONS>;
chomp(@libs); chomp(@exons);
close(LIBS); close(EXONS);

foreach my $lib (@libs) {

    my $reffil = $assemdir . $lib . "/" . $lib . "_refs.fasta";
    open REFS, ">$reffil" or die "cannot the ref fasta file";

    foreach my $exonfile (@exons) {

        if ($exonfile =~ /(ENS\S+)_(exon\d+)_/) {

            my $exon = $1 . "_" . $2;
            my $bestcontigfil = $assemdir . $lib . "/" . $exon . ".fa.best";

            if(-e $bestcontigfil) {

                open CTG, "<$bestcontigfil" or die "could not open the contig file";
                my @contig = <CTG>;
                close(CTG); 
                chomp(@contig);
                my $seq = $contig[1];
                print REFS ">" . $exon ."\n". $seq ."\n"; 

            }
        }
    }
}
