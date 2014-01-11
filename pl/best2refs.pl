use warnings;
use strict;

# Gather all the reads from the .best file for each exon
# into a file with _refs.fasta for each library

# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $libfil   = $assemdir . "libraries.txt";
# my $exonlist = $assemdir . "refs/targets/targetexons.txt.all";

my ($assemdir, $libfil, $exonlist) = @ARGV;

open LIBS, "<$libfil" or die "could not open the lib file";
open EXONS, "<$exonlist" or die "could not open the lib file";

my @libs  = <LIBS>;
my @exons = <EXONS>;
chomp(@libs); chomp(@exons);
close(LIBS); close(EXONS);

my $gathercontigs_dir = "$assemdir/gathercontigs_by_target_exon/";

foreach my $lib (@libs) {

    my $reffil = "$assemdir/$lib/${lib}_best2refs.fasta";
    open REFS, ">$reffil" or die "Could not open the output best2refs fasta file $reffil\n";

    foreach my $exon (@exons) {

        if ($exon =~ /(ENS\S+)_(exon\S+)/) { 

            my $prot = $1; 
            my $exon_name = $1 . "_" . $2;

            my $bestcontig_distrib_dir = "$assemdir/$lib/$prot/${prot}_bestcontig_distrib/";
            my $bestcontigfil = "$bestcontig_distrib_dir/${exon_name}_velvet_contigs.cap3ed.exonerated.filtered.best_contig.fasta";      

            open CTG, "<$bestcontigfil" or die "Could not open the best contig file at $bestcontigfil\n";
            my @contig = <CTG>;
            chomp(@contig);
            close(CTG); 

            if (scalar(@contig) == 2) {
                my $seq = $contig[1];
                print REFS ">" . $exon_name ."\n". $seq ."\n"; 
            } elsif (scalar(@contig) == 0) {
                # This is fine
                # warn "[WARNING best2refs] No best contig for $exon_name\n";
            } else {
                die "Incorrect number of contigs/contig format in $bestcontigfil\n";
            }

        } else {
            die "Invalid exon name $exon found in $exonlist\n";
        }
        
    }
}
