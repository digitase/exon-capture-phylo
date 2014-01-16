use warnings;
use strict;

# Gather all the reads from the .best file for each exon
# into a file with _refs.fasta for each library

# my $assemdir = "/home2/jgb/assemble/crypto/";
# my $libfil   = $assemdir . "libraries.txt";
# my $exonlist = $assemdir . "refs/targets/targetexons.txt.all";

my ($lib, $assemdir, $exonlist) = @ARGV;

my $best2refs_dir = "$assemdir/$lib/${lib}_best2refs/";
unless(-e $best2refs_dir or mkdir $best2refs_dir) { die "Unable to create $best2refs_dir\n"; }

open EXONS, "<$exonlist" or die "could not open the lib file";
my @exons = <EXONS>;
close(EXONS);
chomp(@exons);

my $gathercontigs_dir = "$assemdir/gathercontigs_by_target_exon/";

my $reffil = "$best2refs_dir/${lib}_best2refs.fasta";
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
            # This is fine.
            # warn "[WARNING best2refs] No best contig for $exon_name\n";
        } else {
            die "Incorrect number of contigs/contig format in $bestcontigfil\n";
        }

    } else {
        die "Invalid exon name $exon found in $exonlist\n";
    }
    
}
