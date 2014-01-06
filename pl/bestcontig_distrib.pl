use warnings;
use strict;

# # Carlia exons identified by homology to anolis proteins
# # From a carlia transcriptome, we identify carlia exons and create probes by homology to anolis target proteins
# # these are simple fasta files
# my $exondir  = "/home2/jgb/assemble/crypto/refs/targets/";

# # Format:
# # ENSACAP00000000005_exon1_Carlia.fasta
# my $exonlist = "/home2/jgb/assemble/crypto/refs/targets/targetexons.txt.all";

# # sample names list
# my $assemdir = "/home2/jgb/camaxlibs/";
# my $libfil   = $assemdir . "camaxalllibs.txt";

# # here we have the anolis protein target sequences to blast to
# my $protfil  = "/home2/jgb/assemble/crypto/anolisproteinlist.txt";
# my $protseqs = "/home2/jgb/assemble/crypto/refs/anolistargetproteins/";

# # Blast to anolis proteins instead of carlia exon dnaseq to avoid biasing for exons based on divergence to carlia
# # We blast to all anolis proteins.
# # In the case that the best hit when blasting our assembled by prot exon is not the exon that we assembled by, we ignore it
# # TODO add database creation to script
# my $blastdb  = "/home2/jgb/blastdb/anolis_carolinensis/Anolis_carolinensis.AnoCar2.0.67.pep.all.fa";

# take from argv instead
my ($lib, $assemdir, $libfil, $exondir, $exonlist, $protseqs, $protfil, $all_target_seqs, $blastdb, $minoverlap) = @ARGV;

# Create the target BLAST database unless it already exists
unless(-e "$blastdb.pin") {
    system("makeblastdb -dbtype prot -in $all_target_seqs -out $blastdb");
}

open EXONS, "<$exonlist" or die "could not open the lib file";
my @exons = <EXONS>;
chomp(@exons);
close(EXONS);

foreach my $exonfile (@exons) {

    #sleep(2);

# exonlist entries are of the form ENSACAP00000021611_exon1_Sapro.fasta
    if ($exonfile =~ /(ENS\S+)_(exon\d+)_/) { 

# Grab protein and exon name
        my $prot = $1; 
        my $exon = $1 . "_" . $2;   

# Grab exonerate alignment region bounds
        my ($lower, $upper) = getlimits($exondir . $exonfile);

# Existing file with exonerate header, cap3 contigs and any uncatted velvet contigs
        my $contigsallkexonerate = $assemdir . $lib . "/" . $prot . "_velvetsixk.fa.cap3out.exonerate";
# File with the original velvet contigs
        my $exonlibfil = $assemdir . $lib . "/" . $exon . ".fa";      

# Prefilter
# Clip length > 65%
# Clip to intron-exon boundaries
        my $call1 = parseexon($contigsallkexonerate, $exonlibfil, $lower, $upper);

        my $call2 = performRBH($exonlibfil, $prot, $blastdb);

   } else {
      print "could not recognise the exon...";
   }

} # end foreach @exons


sub getlimits {

    # TODO This assumes the exon has been pre-exonerated
    my $exfilex = $_[0] . ".exonerate";
    open EFE, "<$exfilex" or die "could not open exonfile";
    my @lines = <EFE>;
    my $nameline = $lines[2];

    # exonerate --ryo >%ti b%qab e%qae p%pi\n%tas\n
    # >ENSACAP00000013015_exon4_Carlia b157 e231 p81.08
    # ACACGAGCCACCTCCCGGGAGGAGAAGAACCTTCAAAGCTTTCTGGAGCACCCAAAGGAGAAGTGGGTAG
    # AGAGTGCCTTTGAGGTGGACGGGCCACACTACTATATAGTCATGGCACTCCACATCCTGCCCCCGGAGAG
    # GTGGAAAGCCATGCGCATCGACATCCTCAGGCGGCTGCTGGTGATCTCCCAGGCCCGGGTGGTGTCTCCA
    # GGGGGAGCAAGC

    # Get beginning and end of the query region in the alignment
    if ($nameline =~ / b(\d+) e(\d+) p/) {
        my $b = $1; my $e = $2; return($b, $e);
    } else {
        my $b = "0"; my $e = "0"; return($b, $e);
    }

}


sub parseexon {

    my $exonerated = $_[0];
    my $unexonerated = $_[1];
    my $lower = $_[2];
    my $upper = $_[3];

    open EX,   "<$exonerated" or die "could not open exonfile $exonerated";

    #remove exonerate lines
    my @filtlines = ();
    my @lines = <EX>;

    foreach my $line (@lines) {

    #TODO Surely this is equivalent to using NOT
    # Remove exonerate lines and blank lines
       if ($line =~ /^Command|^Host|completed exonerate|^\n$/) {
           # do nothing
       } else {
           push(@filtlines, $line);
       }

    }
    
    # check overlap with exon and print
    my $newfilstring = join('', @filtlines);
    # Split into records by '>'
    my @tmplines = split(/\n>/, $newfilstring);

    open UNEX, ">$unexonerated" or die "could not open exonfile $unexonerated";
    my $t = 1;
    foreach my $tmplin (@tmplines) {

        my ($name, $tmpseq) = split(/\n/,$tmplin, 2);
        $name =~ s/>//;

        my $b; my $e; 
        if ($name =~ / b(\d+) e(\d+) p/) {
            $b = $1;
            $e = $2;
        } else {
            $b = "0";
            $e = "0";
        }

        # is there any overlap 
        if (($b >= $lower && $b < $upper) || ($e > $lower && $e <= $upper) || ($b <= $lower && $e >= $upper)) {
  
            # is there enough overlap
            my $bminuslower = 0;
            my $upperminuse = 0;

            if ($b > $lower) { $bminuslower = $b-$lower };
            if ($upper > $e) { $upperminuse = $upper-$e };

            my $overlap =  ($upper - $lower - $bminuslower - $upperminuse  + 1) / ($upper - $lower + 1);

            if ( $overlap > $minoverlap ){
      
                $tmpseq =~ s/\n//g;
                my $newname = "excontig$t";
                $t++;
                print UNEX ">$newname\n$tmpseq\n";
            }
        }
    }
}


sub performRBH {
# Blast the contig against all anolis proteisn
    my ($exonlibfilclust, $prot, $blastdb) = @_;
    my $blastout = "$exonlibfilclust.blast";
    my $bestout = "$exonlibfilclust.best";
    system("blastall -i $exonlibfilclust -p blastx -d $blastdb -o $blastout -m 8 -e 1E-10");

    open BLAST, "<$blastout" or die "could not open exonfile";
    my @blastlines = <BLAST>;
    chomp(@blastlines);
    my $maxbitscore = 0;
    my $bestprothit = "";
    my $bestcontig  = "";

    foreach my $line (@blastlines) {
       my @linbits = split(/\t/, $line);
       if ($linbits[11] > $maxbitscore) {
               $maxbitscore = $linbits[11];
               $bestprothit = $linbits[1];
               $bestcontig  = $linbits[0];
       }
    }

# Make sure the best hit is indeed the exon the contig was assembled off    
    if ($bestprothit eq $prot){

       #print "$bestprothit\t$prot\t$bestcontig\n" ;
       open BEST,    ">$bestout" or die "could not open exonfile";
       open CONTIGS, "<$exonlibfilclust" or die "could not open exonfile";
       my @contiglines = <CONTIGS>;
       my $contigstring = join('',@contiglines);
       my @tmplines = split(/\n>/,$contigstring);

       foreach my $tmplin (@tmplines) {

          my ($name, $tmpseq) = split(/\n/,$tmplin, 2);
          $name =~ s/>//;

# Then we can write out the best contig to the best file
          if ($name eq $bestcontig) {
              $tmpseq =~ s/\n//g;
              print BEST ">$name\n$tmpseq\n";
          }

       }

    }

}
  
