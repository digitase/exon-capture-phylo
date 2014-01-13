use warnings;
use strict;

use List::Util qw(max min);

my ($lib, $assemdir, $all_exons_file, $exonlist, $blastdb, $minoverlap) = @ARGV;

# Create the target BLAST database unless it already exists
my $blast_dbs_dir = $assemdir . "/blast_dbs/";
unless(-e "$blast_dbs_dir/$blastdb.pin") {
    system("makeblastdb -dbtype prot -in $assemdir/all_proteins.fasta -out $blast_dbs_dir/$blastdb");
}

# Grab target exon IDs
open EXONS, "<$exonlist" or die "Could not open list of target exons file $exonlist";
my @exons = <EXONS>;
chomp(@exons);
close(EXONS);

my $assemlib = "$assemdir/$lib/";

foreach my $exon (@exons) {

    # exonlist entries are of the form ENSACAP00000021611_exon1_Sapro
    # TODO change to alphanumberic regex
    if ($exon =~ /(ENS\S+)_(exon\S+)/) { 

        my $prot = $1; 
        my $exon_name = $1 . "_" . $2;   

        my $bestcontig_distrib_dir = "$assemlib/$prot/${prot}_bestcontig_distrib/";
        unless(-e $bestcontig_distrib_dir or mkdir $bestcontig_distrib_dir) { die "Could not make $bestcontig_distrib_dir\n"; }

        # Get exon sequence from all exons file
        my $exonerate_query = "$assemlib/$prot/${prot}_catcontigs/$prot.fasta";
        my $exonerate_target = "$bestcontig_distrib_dir/$exon_name.fasta";
        system("perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($exon) } print if \$c' $all_exons_file > $exonerate_target");

        my $exonerate_out = "$bestcontig_distrib_dir/$exon_name.exonerated.fasta";
        # get overlap between target exon and anolis
        my ($region_start, $region_end) = getTargetRegionInProtein($exonerate_query, $exonerate_target, $exonerate_out);

        # get overlaps between assembled contigs and anolis, then filter for those contigs with good overlap with target region
        # File from catcontigs
        my $exonerated_contigs = "$assemlib/$prot/${prot}_catcontigs/${prot}_velvet_contigs.cap3ed.exonerated.fasta";
        # output file with contigs that overlap to the current exon
        my $filtered_contigs = "$bestcontig_distrib_dir/${exon_name}_velvet_contigs.cap3ed.exonerated.filtered.fasta";      
        filterExoneratedContigs($exonerated_contigs, $filtered_contigs, $region_start, $region_end, $exon_name, $minoverlap);

        # determine best contig for the exon by blastx
        getBestContig($filtered_contigs, $prot, $blastdb, $blast_dbs_dir, $exon_name);

    } else {
        print "Exon naming format incorrect for: $exon\n";
    }

    # sleep(2);
}

sub getTargetRegionInProtein {
    my ($exonerate_query, $exonerate_target, $exonerate_out) = @_;

    my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
    system("exonerate --model protein2genome --query $exonerate_query --target $exonerate_target --showvulgar no --showalignment no --ryo $ryo > $exonerate_out");
    
    # Remove exonerate lines
    open EFE, "<$exonerate_out" or die "Could not open exonerate output file $exonerate_out\n";
    my @lines = grep(!/^Command line:|^Hostname:|-- completed exonerate analysis|^\n$/, <EFE>);
    close EFE;

    # Get beginning and end of the query region in the alignment
    my $nameline = $lines[0];
    if (scalar(@lines) == 0) {
        warn "[WARNING bestcontig_distrib] No prot-exon exonerate alignment between $exonerate_query and $exonerate_target\n";
        return("0", "0");
    } elsif ($nameline =~ / b(\d+) e(\d+) p/) {
        my $b = $1; my $e = $2;
        return($b, $e);
    } else {
        die "Invalid exonerate sequence ID line $nameline in $exonerate_out\n";
    }
}

# Check for assembled contigs to overlap with exon
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
        if (scalar(@contig_file_lines) == 0) {
            warn "[WARNING bestcontig_distrib] No prot-contig exonerate alignment in $exonerated_contigs detected\n";
            $contig_name = ""; $b = "0"; $e = "0";
        } elsif ($contig_name_line =~ /^>(\S+) b(\d+) e(\d+) p/) {
            $contig_name = $1; $b = $2; $e = $3;
        } else {
            die "Invalid exonerate sequence ID line $contig_name_line in $exonerated_contigs\n";
        }

        # Move on if there is no overlap
        my $overlap = max (0, ((min ($region_end, $e)) - (max ($region_start, $b))));
        next unless $overlap;

        # Check overlap proportion
        my $overlap_ratio = $overlap / ($region_end - $region_start);
        if ($overlap_ratio >= $minoverlap) {
            $contig_seq =~ s/\n//g;
            print OUT ">${exon_name}_contig_$contig_num\n$contig_seq\n";
            $contig_num++;
        } else {
            warn "[WARNING bestcontig_distrib] $contig_name failed filtering. Overlap=$overlap_ratio. Required:$minoverlap\n";
        }
    }
    close OUT;
}

sub getBestContig {
    # Blast the contig against all anolis proteisn
    my ($filtered_contigs_file, $prot, $blastdb, $blast_dbs_dir, $exon_name) = @_;

    my $filtered_contigs_basename = $filtered_contigs_file;
    $filtered_contigs_basename =~ s/\.fasta//;

    my $blastout = "$filtered_contigs_basename.against_all.blast";
    my $bestout = "$filtered_contigs_basename.best_contig.fasta";
    system("blastall -i $filtered_contigs_file -p blastx -d $blast_dbs_dir/$blastdb -o $blastout -m 8 -e 1E-10");
    
    open BLAST, "<$blastout" or die "Could not open blastx output file $blastout\n";
    my @blastlines = <BLAST>;
    chomp(@blastlines);
    my $maxbitscore = 0;
    my $bestprothit = "";
    my $bestcontig  = "";

    # Get info for the best blast hit score
    foreach my $line (@blastlines) {
       my @linbits = split(/\t/, $line);
       if ($linbits[11] > $maxbitscore) {
               $bestcontig  = $linbits[0];
               $bestprothit = $linbits[1];
               $maxbitscore = $linbits[11];
       }
    }

    # Make sure the best hit is indeed the exon the contig was assembled off    
    open BEST, ">$bestout" or die "Could not open best contig output file $bestout\n";
    if ($bestprothit eq $prot) {
        open CONTIGS, "<$filtered_contigs_file" or die "Could not open filtered assembled contigs file $filtered_contigs_file\n";
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
    } else {
        warn "[WARNING bestcontig_distrib] The best contig $bestcontig failed reciprocal best-hit blast. Best hit was to $bestprothit\n";
    }
}
  
