use warnings;
use strict;

use List::Util qw(max min);

my ($lib, $assemdir, $libfil, $exondir, $exonlist, $all_target_seqs, $blastdb, $minoverlap) = @ARGV;

# Create the target BLAST database unless it already exists
my $blast_dbs_dir = $assemdir . "/blast_dbs/";
chdir("$blast_dbs_dir") or die "Cannot chdir to $blast_dbs_dir\n";
unless(-e "$blastdb.pin") {
    system("makeblastdb -dbtype prot -in $all_target_seqs -out $blastdb");
}
chdir("$assemdir") or die "Cannot chdir to $assemdir\n";

# Grab target exon IDs
open EXONS, "<$exonlist" or die "could not open the lib file";
my @exons = <EXONS>;
chomp(@exons);
close(EXONS);

my $assemlib = "$assemdir/$lib/";

foreach my $exonfile (@exons) {

    # exonlist entries are of the form ENSACAP00000021611_exon1_Sapro.fasta
    # TODO change to alphanumberic regex
    if ($exonfile =~ /(ENS\S+)_(exon\d+)_/) { 

        my $prot = $1; 
        my $exon_name = $1 . "_" . $2;   

        my $bestcontig_distrib_dir = "$assemlib/$prot/${prot}_bestcontig_distrib/";
        unless(-e $bestcontig_distrib_dir or mkdir $bestcontig_distrib_dir) { die "Could not make $bestcontig_distrib_dir\n"; }

        # TODO add this to the pipeline; combine all exons into a file, don't use a dir
        # get overlap between target exon and anolis
        my ($region_start, $region_end) = getTargetRegionInProtein($exondir . $exonfile);

        # get overlaps between assembled contigs and anolis, then filter for those contigs with good overlap with target region
        # File from catcontigs
        my $exonerated_contigs = "$assemlib/$prot/${prot}_catcontigs/${prot}_velvet_contigs.cap3ed.exonerated.fasta";
        # output file
        my $filtered_contigs = "$bestcontig_distrib_dir/$exon_name.good_overlap.fasta";      
        filterExoneratedContigs($exonerated_contigs, $filtered_contigs, $region_start, $region_end, $exon_name, $minoverlap);

        # determine best contig for the exon by blastx
        getBestContig($filtered_contigs, $prot, $blastdb, $blast_dbs_dir, $exon_name, $minoverlap);

    } else {
        print "Exon naming format incorrect for: $exonfile\n";
    }

    # sleep(2);
}

sub getTargetRegionInProtein {
    # TODO This assumes the exon has been pre-exonerated
    my $exfilex = $_[0] . ".exonerate";
    
    open EFE, "<$exfilex" or die "could not open exonfile";
    my @lines = <EFE>;
    close EFE;
    my $nameline = $lines[2];

    # Get beginning and end of the query region in the alignment
    if ($nameline =~ / b(\d+) e(\d+) p/) {
        my $b = $1; my $e = $2;
        return($b, $e);
    } else {
        die "Invalid exonerate sequence ID line $nameline in $exfilex\n";
    }
}


sub filterExoneratedContigs {
    my ($exonerated_contigs, $filtered_contigs, $region_start, $region_end, $exon_name, $minoverlap) = @_;

    # Remove exonerate lines
    open IN, "<$exonerated_contigs" or die "Could not open exonerated contigs file $exonerated_contigs\n";
    my @contig_file_lines = grep(!/^Command line:|^Hostname:|-- completed exonerate analysis|^\n$/, <IN>);
    close IN;

    # Split into records by '>'
    my @contigs = split(/\n(?=>)/, join("", @contig_file_lines));

    my $contig_num = 0;
    open OUT, ">$filtered_contigs" or die "Could not open filtered contigs output file $filtered_contigs\n";
    foreach my $contig (@contigs) {
        my ($contig_name_line, $contig_seq) = split(/\n/, $contig, 2);

        my $contig_name; my $b; my $e; 
        if ($contig_name_line =~ /^>(\S+) b(\d+) e(\d+) p/) {
            $contig_name = $1; $b = $2; $e = $3;
        } else {
            die "Invalid exonerate sequence ID line $contig_name_line in $exonerated_contigs\n";
        }

        my $overlap = max (0, ((min ($region_end, $e)) - (max ($region_start, $b))));
        my $overlap_ratio = $overlap / ($region_end - $region_start);

        if ($overlap >= $minoverlap) {
            $contig_seq =~ s/\n//g;
            print OUT ">${exon_name}_contig_$contig_num\n$contig_seq\n";
            $contig_num++;
        }
    }
    close OUT;
}


sub getBestContig {
    # Blast the contig against all anolis proteisn
    my ($filtered_contigsclust, $prot, $blastdb, $blast_dbs_dir) = @_;
    my $blastout = "$filtered_contigsclust.against_all.blast";
    my $bestout = "$filtered_contigsclust.best_contig";
    system("blastall -i $filtered_contigsclust -p blastx -d $blast_dbs_dir/$blastdb -o $blastout -m 8 -e 1E-10");
    
    open BLAST, "<$blastout" or die "could not open exonfile";
    my @blastlines = <BLAST>;
    chomp(@blastlines);
    my $maxbitscore = 0;
    my $bestprothit = "";
    my $bestcontig  = "";

    foreach my $line (@blastlines) {
       my @linbits = split(/\t/, $line);
       if ($linbits[11] > $maxbitscore) {
               $bestcontig  = $linbits[0];
               $bestprothit = $linbits[1];
               $maxbitscore = $linbits[11];
       }
    }

    # Make sure the best hit is indeed the exon the contig was assembled off    
    open BEST,    ">$bestout" or die "could not open exonfile";
    if ($bestprothit eq $prot) {
        open CONTIGS, "<$filtered_contigsclust" or die "could not open exonfile";
        my @tmplines = split(/\n>/, join('', <CONTIGS>));
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
  
