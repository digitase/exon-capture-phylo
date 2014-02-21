use strict;
use warnings;

# Sample name, samples directory, output directory, FASTA with all the targets to make the blast db, BLAST database name
my ($lib, $readdir, $assemdir, $fwd_suffix, $rev_suffix, $unpaired_suffix, $target_seqs_list, $np, $blastall_path) = @ARGV;

# Expectation value for blastx
my $eval = "1e-9";

my $blast_dbs_dir = $assemdir . "/blast_dbs/";
filtAssemb($readdir, $assemdir, $blast_dbs_dir, $lib, $fwd_suffix, $rev_suffix, $unpaired_suffix, $eval, $np);

sub filtAssemb {
    my ($readdir, $assemdir, $blast_dbs_dir, $lib, $fwd_suffix, $rev_suffix, $unpaired_suffix, $eval, $np) = @_;

    # Create directory for sample
    my $assemlib = $assemdir . $lib . "/";
    unless(-d $assemlib or mkdir $assemlib) { die "Could not mkdir $assemlib\n"; }    

    # Paths to read files
    # This section determines what the suffixes of the input read files need to be
    # TODO consider generic suffixes    

    my $fil_1 = "$readdir/$lib" . "_$fwd_suffix.fastq.gz";
    my $fil_2 = "$readdir/$lib" . "_$rev_suffix.fastq.gz";
    my $fil_u = "$readdir/$lib" . "_$unpaired_suffix.fastq.gz";

    my $blast_1 = "$assemlib/$lib" . "_$fwd_suffix.against_targets.blast";
    my $blast_2 = "$assemlib/$lib" . "_$rev_suffix.against_targets.blast";
    my $blast_u = "$assemlib/$lib" . "_$unpaired_suffix.against_targets.blast";

    my $fasta_1 = "$assemlib/$lib" . "_$fwd_suffix.fasta";
    my $fasta_2 = "$assemlib/$lib" . "_$rev_suffix.fasta";
    my $fasta_u = "$assemlib/$lib" . "_$unpaired_suffix.fasta";

    # do the blasting
    unless(-e "$blast_1") { blastProts($fil_1, $blast_1, $fasta_1, $blast_dbs_dir, $eval, $np); }
    unless(-e "$blast_2") { blastProts($fil_2, $blast_2, $fasta_2, $blast_dbs_dir, $eval, $np); }
    unless(-e "$blast_u") { blastProts($fil_u, $blast_u, $fasta_u, $blast_dbs_dir, $eval, $np); }

    # for each exon that was hit by reads from a file, collate those reads
    my $call_1  = getbest($assemlib, $fasta_1, $blast_1, "1" , $target_seqs_list);
    my $call_2  = getbest($assemlib, $fasta_2, $blast_2, "2" , $target_seqs_list);
    my $call_u  = getbest($assemlib, $fasta_u, $blast_u, "u" , $target_seqs_list);

    # for each exon that was hit by reads from a file, collate those reads with the same sequence ID from the paired file
    my $call_1p = getbest($assemlib, $fasta_2, $blast_1, "1p", $target_seqs_list);
    my $call_2p = getbest($assemlib, $fasta_1, $blast_2, "2p", $target_seqs_list);
}


sub blastProts {
    my ($fil, $blast, $fasta, $blast_dbs_dir, $eval, $np) = @_;
    
    my $blast_call = "$blastall_path -p blastx -d $blast_dbs_dir/target_proteins -e $eval -m 8 -I T";
    my $start_time = localtime();
    print "Blasting $fil against target proteins with legacy BLAST at $start_time\n";

    system(qq(
        zcat $fil | 
        awk '{if(NR % 4 == 1 || NR % 4 == 2) {sub(/@/, ">"); print; } }' | 
        tee $fasta |
        parallel -j $np --block 200K --recstart '>' --pipe $blast_call > $blast
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

=head1 NAME

assembleByProt

=head1 SYNOPSIS

BLASTx sample reads against target proteins and gather reads aligned to each target

=head1 USAGE

$lib, $readdir, $assemdir, $target_seqs_list, $np
Sample name, samples directory, output directory, database name, target protein IDs list, num of parallel blast processes

=head1 DESCRIPTION OF ARGUMENTS

to be continued...

=head1 SUBROUTINES

Use CPAN style http://juerd.nl/site.plp/perlpodtut

=head1 DIAGNOSTICS

List of errors...

=head1 DEPENDENCIES

NCBI blastall
GNU parallel

=head1 BUGS AND LIMITATIONS

NCBI blastall is being phased out...

=head1 AUTHORS

Dr. Jason Bragg, ANU
Ben Bai, ANU

=head1 AUTHORS' COMMENTS

=head2 Why not BLAT?

BLAT (BLAST-like Alignment Tool) is a sequence alignment tool similar to BLAST but structured differently. BLAT quickly finds similarity in DNA and protein but it needs an exact or nearly-exact match to find a hit. Therefore Blat is not as flexible as BLAST. Since BLAST can find much more remote matches than Blat, it is the recommended tool when searching more distantly related sequences.

=head2 Why not BLAST+?

Speed...

=head2 Why the stringent expectation value?

Blah...

=head1 LICENCE AND COPYRIGHT

GNU Public

=cut
