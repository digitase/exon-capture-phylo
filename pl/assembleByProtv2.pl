# blastx reads against targets
# reads with significant hits to each exon per file
# BLAT (BLAST-like Alignment Tool) is a sequence alignment tool similar to BLAST but structured differently. BLAT quickly finds similarity in DNA and protein but it needs an exact or nearly-exact match to find a hit. Therefore Blat is not as flexible as BLAST. Since BLAST can find much more remote matches than Blat, it is the recommended tool when searching more distantly related sequences.

use strict;
use warnings;

# Sample name, samples directory, output directory, FASTA with all the targets to make the blast db, BLAST database name
my ($lib, $readdir, $assemdir, $all_prot_seqs, $adb, $target_seqs_list) = @ARGV;

# # Copy in all prots file, using the first field as the ID
# system("awk '{print \$1}' $all_prot_seqs > $assemdir/all_proteins.fasta");

# # Create blast database directory
# unless(-d $blast_dbs_dir or mkdir $blast_dbs_dir) {
    # die "Could not create blast database output directory $blast_dbs_dir\n";
# }

# # Extract sequences in target list from proteins file
# my $target_seqs = "$assemdir/target_proteins.fasta";
# system("perl -ne 'if (/^>(\\S+)/) {\$c=\$i{\$1}}\$c?print:chomp;\$i{\$_}=1 if \@ARGV' $target_seqs_list $assemdir/all_proteins.fasta > $target_seqs");

# # Create the target BLAST database unless it already exists
# unless(-e "$blast_dbs_dir/$adb.pin") {
    # system("makeblastdb -dbtype prot -in $target_seqs -out $blast_dbs_dir/$adb");
# }

# Expectation value for blastx
my $eval = "1e-9";
# Number of parallel blast processes to use
my $np = "16";
# Use blastall blastx instead of blast+ blastx?
my $use_legacy_blast = 1;

my $blast_dbs_dir = $assemdir . "/blast_dbs/";
filtAssemb($readdir, $assemdir, $blast_dbs_dir, $lib, $adb, $eval, $np, $use_legacy_blast);

sub filtAssemb {
    my ($readdir, $assemdir, $blast_dbs_dir, $lib, $adb, $eval, $np, $use_legacy_blast) = @_;

    # Create directory for sample
    my $assemlib = $assemdir . $lib . "/";
    unless(-d $assemlib or mkdir $assemlib) { die "Could not mkdir $assemlib\n"; }    

    # Paths to read files
    # This section determines what the suffixes of the input read files need to be
    # TODO consider generic suffixes    

    my $fil_1 = "$readdir/$lib" . "_1_final.fastq.gz";
    my $fil_2 = "$readdir/$lib" . "_2_final.fastq.gz";
    my $fil_u = "$readdir/$lib" . "_u_final.fastq.gz";

    my $blast_1 = "$assemlib/$lib" . "_1_final.against_targets.blast";
    my $blast_2 = "$assemlib/$lib" . "_2_final.against_targets.blast";
    my $blast_u = "$assemlib/$lib" . "_u_final.against_targets.blast";

    my $fasta_1 = "$assemlib/$lib" . "_1_final.fasta";
    my $fasta_2 = "$assemlib/$lib" . "_2_final.fasta";
    my $fasta_u = "$assemlib/$lib" . "_u_final.fasta";

    # do the blasting
    unless(-e "$blast_1") { blastProts($fil_1, $blast_1, $fasta_1, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }
    unless(-e "$blast_2") { blastProts($fil_2, $blast_2, $fasta_2, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }
    unless(-e "$blast_u") { blastProts($fil_u, $blast_u, $fasta_u, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }

    # for each exon that was hit by reads from a file, collate those reads
    my $call_1  = getbest($assemlib, $fasta_1, $blast_1, "1" , $target_seqs_list);
    my $call_2  = getbest($assemlib, $fasta_2, $blast_2, "2" , $target_seqs_list);
    my $call_u  = getbest($assemlib, $fasta_u, $blast_u, "u" , $target_seqs_list);

    # for each exon that was hit by reads from a file, collate those reads with the same sequence ID from the paired file
    my $call_1p = getbest($assemlib, $fasta_2, $blast_1, "1p", $target_seqs_list);
    my $call_2p = getbest($assemlib, $fasta_1, $blast_2, "2p", $target_seqs_list);
}

sub blastProts {
    my ($fil, $blast, $fasta, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast) = @_;
    
    my $blast_call;
    my $start_time = localtime();

    if($use_legacy_blast) {
        print "Blasting $fil against $adb with legacy BLAST at $start_time\n";
        $blast_call = "blastall -p blastx -d $blast_dbs_dir/$adb -e $eval -m 8 -I T"
    } else {
        print "Blasting $fil against $adb with BLAST+ at $start_time\n";
        $blast_call = "blastx -db $blast_dbs_dir/$adb -evalue $eval -outfmt 6 -show_gis"
    }

    system(qq(
        zcat $fil | 
        awk '{if(NR % 4 == 1 || NR % 4 == 2) {sub(/@/, ">"); print; } }' | 
        tee $fasta |
        parallel -j $np --block 1M --recstart '>' --pipe $blast_call > $blast
    ));
}

#  
sub getbest {
    my ($assemlib, $fasta, $blast, $f, $target_seqs_list) = @_;

    open TARGET_SEQS_LIST, "<$target_seqs_list" or die "could not open the target IDs list";
    my @protnames = <TARGET_SEQS_LIST>;
    chomp(@protnames);
    close(TARGET_SEQS_LIST);

    # Hash of sequence names -> sequences
    my %fasta;
    open(FA, "<$fasta") or die "Failed to open FASTA file $fasta\n";
    while (<FA>)
    {
        my $seqnm = $_;
        # This is fine as read files have alternating ID read lines
        my $seq   = <FA>;
        $seqnm =~ s/>//;
        chomp($seqnm); chomp($seq);
        $fasta{ $seqnm } = $seq;
    }
    close FA;
        
    # Hash of protein names -> hit sequence names to the prot -> sequences
    my %prothits = map { $_ => {} } @protnames;
    open(BLOUT, "<$blast") or die "Failed to open BLAST output file $blast\n";
    while(<BLOUT>) 
    {
       my @linbits = split(/\t/, $_);
       my $hitseqnm = $linbits[0];
       my $hitsprot = $linbits[1];
       $prothits{$hitsprot}{$hitseqnm} = $fasta{ $hitseqnm };
    }

    foreach my $prot (keys %prothits) {

        # Create script output directory structure for target
        my $poutdir = "$assemlib/$prot/";
        unless(-d $poutdir or mkdir $poutdir) { die "Could not create exon output directory $poutdir\n"; }

        my $assemble_by_prot_dir = "$poutdir/${prot}_assemble_by_prot/";
        unless(-d $assemble_by_prot_dir or mkdir $assemble_by_prot_dir) {
            die "Could not create assembleByProt output directory $assemble_by_prot_dir\n";
        }

        # Gather reads that blast hit onto target
        my $poutfil = "$assemble_by_prot_dir/${prot}_${f}_hitreads.fasta";
        open(POUT, ">$poutfil") or die "Failed to open poutfil $poutfil\n"; 
        foreach my $seqhitprot (keys %{$prothits{$prot}}) {
            print POUT ">$seqhitprot\n$prothits{$prot}{$seqhitprot}\n";
        }
        close(POUT);
    }
}


__END__


