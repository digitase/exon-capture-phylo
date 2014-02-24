#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

use List::Util qw(max min);

my ($lib, $assemdir, $all_exons_file, $exonlist, $minoverlap, $exonerate_path, $blastall_path) = @ARGV;

my $assemlib = "$assemdir/$lib/";
my $blast_dbs_dir = "$assemdir/blast_dbs/";

# Grab target exon IDs
open EXONS, "<$exonlist" or die "[ERROR bestcontig_distrib_dir $lib] Could not open list of target exons file $exonlist";
my @exons = <EXONS>;
chomp(@exons);
close(EXONS);

foreach my $exon_name (@exons) {
    if ($exon_name =~ /^(\S+?)_\S+/) {

        # Non-whitespace portion of exon ID before underscore should be the orthologous protein ID
        my $prot = $1; 

        my $bestcontig_distrib_dir = "$assemlib/$prot/${prot}_bestcontig_distrib/";
        unless(-e $bestcontig_distrib_dir or mkdir $bestcontig_distrib_dir) { 
            die "[ERROR bestcontig_distrib_dir $lib] Could not make $bestcontig_distrib_dir\n";
        }

        # Get exon sequence from all exons file
        my $exonerate_query = "$assemlib/$prot/${prot}_catcontigs/$prot.fasta";
        my $exonerate_target = "$bestcontig_distrib_dir/$exon_name.fasta";
        system("perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($exon_name) } print if \$c' $all_exons_file > $exonerate_target");

        # Get overlap region bounds between target exon and target protein with exonerate
        my $exonerate_out = "$bestcontig_distrib_dir/$exon_name.exonerated.fasta";
        my ($region_start, $region_end) = getTargetRegionInProtein($exonerate_query, $exonerate_target, $exonerate_out);

        # Read overlap regions for each assembled contig and compare to the above region for length filtering
        my $exonerated_contigs = "$assemlib/$prot/${prot}_catcontigs/${prot}_velvet_contigs.cap3ed.exonerated.fasta";
        my $filtered_contigs = "$bestcontig_distrib_dir/${exon_name}_velvet_contigs.cap3ed.exonerated.filtered.fasta";      
        filterExoneratedContigs($exonerated_contigs, $filtered_contigs, $region_start, $region_end, $exon_name, $minoverlap);

        # Determine best contig for the exon by BLASTx
        getBestContig($filtered_contigs, $prot, $blast_dbs_dir, $exon_name);

    } else {
        warn "[WARNING bestcontig_distrib_dir $lib] Exon naming format incorrect for: $exon_name\n";
    }

    # Uncomment to reduce IO demands
    # sleep(2);
}

# Subroutines

sub getTargetRegionInProtein {
    my ($exonerate_query, $exonerate_target, $exonerate_out) = @_;

    # Find alignment region between target exon and protein
    my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
    system("$exonerate_path --model protein2genome --query $exonerate_query --target $exonerate_target --showvulgar no --showalignment no --ryo $ryo > $exonerate_out");
    
    # Remove exonerate-specific output lines to convert output into FASTA format
    open EFE, "<$exonerate_out" or die "[ERROR bestcontig_distrib_dir $lib] Could not open exonerate output file $exonerate_out\n";
    my @lines = grep(!/^Command line:|^Hostname:|-- completed exonerate analysis|^\n$/, <EFE>);
    close EFE;

    # Return beginning and end of the query region in the alignment by parsing the ryo format
    my $nameline = $lines[0];
    if (scalar(@lines) == 0) {
        warn "[WARNING bestcontig_distrib $lib] No prot-exon exonerate alignment between $exonerate_query and $exonerate_target\n";
        return("0", "0");
    } elsif ($nameline =~ / b(\d+) e(\d+) p/) {
        my $b = $1; my $e = $2;
        return($b, $e);
    } else {
        die "[ERROR bestcontig_distrib_dir $lib] Invalid exonerate sequence ID line $nameline in $exonerate_out\n";
    }
}

sub filterExoneratedContigs {
    my ($exonerated_contigs, $filtered_contigs, $region_start, $region_end, $exon_name, $minoverlap) = @_;

    # Remove exonerate-specific lines from exonerated assembled contigs to convert to FASTA format
    open IN, "<$exonerated_contigs" or die "[ERROR bestcontig_distrib_dir $lib] Could not open exonerated contigs file $exonerated_contigs\n";
    my @contig_file_lines = grep(!/^Command line:|^Hostname:|-- completed exonerate analysis|^\n$/, <IN>);
    close IN;

    # Split FASTA into records
    my @contigs = split(/\n(?=>)/, join("", @contig_file_lines));

    my $contig_num = 0;
    open OUT, ">$filtered_contigs" or die "[ERROR bestcontig_distrib_dir $lib] Could not open filtered contigs output file $filtered_contigs\n";

    # Check each contig
    foreach my $contig (@contigs) {
        my ($contig_name_line, $contig_seq) = split(/\n/, $contig, 2);

        my $contig_name; my $b; my $e; 
        if (scalar(@contig_file_lines) == 0) {
            warn "[WARNING bestcontig_distrib $lib] No prot-contig exonerate alignment in $exonerated_contigs detected\n";
            $contig_name = ""; $b = "0"; $e = "0";
        } elsif ($contig_name_line =~ /^>(\S+) b(\d+) e(\d+) p/) {
            # by extracting contig-protein alignment region
            $contig_name = $1; $b = $2; $e = $3;
        } else {
            die "[ERROR bestcontig_distrib_dir $lib] Invalid exonerate sequence ID line $contig_name_line in $exonerated_contigs\n";
        }

        # and calculating the overlap between that region and the exon-protein region determined in getTargetRegionInProtein
        my $overlap = max (0, ((min ($region_end, $e)) - (max ($region_start, $b))));
        my $overlap_ratio = $overlap / ($region_end - $region_start);

        # and discarding those contigs that do not meet a threshold overlap.
        if ($overlap_ratio >= $minoverlap) {
            $contig_seq =~ s/\n//g;
            print OUT ">${exon_name}_contig_$contig_num\n$contig_seq\n";
            $contig_num++;
        } else {
            # Warn of failures
            warn "[WARNING bestcontig_distrib $lib] $exon_name $contig_name_line failed filtering. Overlap=$overlap_ratio. Required=$minoverlap\n";
        }
    }
    close OUT;
}

sub getBestContig {
    my ($filtered_contigs_file, $prot, $blast_dbs_dir, $exon_name) = @_;

    my $filtered_contigs_basename = $filtered_contigs_file;
    $filtered_contigs_basename =~ s/\.fasta//;

    # BLASTx filtered contigs against ALL reference proteins
    my $blastout = "$filtered_contigs_basename.against_all.blast";
    my $bestout = "$filtered_contigs_basename.best_contig.fasta";
    system("$blastall_path -i $filtered_contigs_file -p blastx -d $blast_dbs_dir/all_proteins -o $blastout -m 8 -e 1E-10");
    
    # Extract output lines
    open BLAST, "<$blastout" or die "[ERROR bestcontig_distrib_dir $lib] Could not open blastx output file $blastout\n";
    my @blastlines = <BLAST>;
    close BLAST;
    chomp(@blastlines);

    my $maxbitscore = 0;
    my $bestprothit = "";
    my $bestcontig  = "";
    
    # Determine alignment with highest bit-score
    foreach my $line (@blastlines) {
       my @linbits = split(/\t/, $line);
       if ($linbits[11] > $maxbitscore) {
               $bestcontig  = $linbits[0]; 
               $bestprothit = $linbits[1];
               $maxbitscore = $linbits[11];
       }
    }

    # Make sure the best hit aligned to the correct target protein
    open BEST, ">$bestout" or die "[ERROR bestcontig_distrib_dir $lib] Could not open best contig output file $bestout\n";
    if ($bestprothit eq $prot) {
        open CONTIGS, "<$filtered_contigs_file" or die "[ERROR bestcontig_distrib_dir $lib] Could not open filtered assembled contigs file $filtered_contigs_file\n";
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
    } elsif (scalar(@blastlines) == 0) {
        warn "[WARNING bestcontig_distrib $lib] No blast hits in $blastout\n";
    } else {
        # Otherwise warn of rejected best hits
        warn "[WARNING bestcontig_distrib $lib] The best contig $bestcontig failed reciprocal best-hit blast. Best hit was to $bestprothit\n";
    }
    close BEST;
}

__END__

=head1 NAME

bestcontig_distrib_dir - Filter merged contigs by comparing assembled length to target exon length by aligning both to the target protein with exonerate and comparing overlap. Then BLASTx all filtered contigs to the reference proteome, and determine representative assemblies for each target exon by selecting the best hit.

=head1 USAGE

=over

=item B<perl bestcontig_distrib.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=item $all_exons_file

Nucleotide FASTA file with all target exon sequences.

=item $exonlist

Text file with target exon IDs.

=item $minoverlap

Length filtering threshold: contig-protein/exon-protein alignment range ratio must >= this value for an assembly to pass.

=item $exonerate_path

Path to exonerate binary.

=item $blastall_path

Path to blastall binary.

=back

=head1 SUBROUTINES

=over

=item getTargetRegionInProtein

Align target protein to target exon and determine alignment beginning and end indices on the protein sequence.

=item filterExoneratedContigs

Align target protein to each assembly and compare alignment range to the return values of getTargetRegionInProtein. Assemblies that fail to meet $minoverlap are discarded.

=item getBestContig

BLASTx filtered assemblies to all reference proteins; select representative sequences for target exons.

=back

=head1 DIAGNOSTICS

=over 

=item [ERROR bestcontig_distrib_dir $lib] Could not make ...

File or directory creation failed. Check that you have adequate permissions in the output directory.

=item [ERROR bestcontig_distrib_dir $lib] Could not open ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=item [WARNING bestcontig_distrib $lib] No prot-exon exonerate alignment ...

A target protein and target exon did not align with exonerate. Check that the exon is named correctly, and is orthologous to the target protein.

=item [ERROR bestcontig_distrib_dir $lib] Invalid exonerate sequence ID line ...

Could not find alignment range in exonerate file due to malformed sequence ID line. Check for file corruption.

=item [WARNING bestcontig_distrib $lib] No prot-contig exonerate alignment ...

Target protein did not align with exonerate to any assembled contigs. Check that previous pipeline stages completed successfully, and that there are assembled contigs in the given filename.

=back

=head1 DEPENDENCIES

exonerate

blastall

=head1 KNOWN BUGS

None.

=head1 NOTES

Exonerate alignment range is reported as the start and end of the alignment region in the query sequence (the target protein). For more information on exonerate "ryo" formats, see https://www.ebi.ac.uk/~guy/exonerate/advanced.html and http://csurs7.csr.uky.edu/cgi-bin/man/man2html?1+exonerate

The regex used to parse target exon names restricts the sequence ID format:
    if ($exon_name =~ /^(\S+?)_\S+/) 

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
