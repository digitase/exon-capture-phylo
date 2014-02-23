# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use strict;
use warnings;

my ($lib, $readdir, $assemdir, $fwd_suffix, $rev_suffix, $unpaired_suffix, $target_seqs_list, $np, $blastall_path, $eval) = @ARGV;

my $blast_dbs_dir = "$assemdir/blast_dbs/";
filtAssemb($readdir, $assemdir, $blast_dbs_dir, $lib, $fwd_suffix, $rev_suffix, $unpaired_suffix, $eval, $np);

# Subroutines

sub filtAssemb {
    my ($readdir, $assemdir, $blast_dbs_dir, $lib, $fwd_suffix, $rev_suffix, $unpaired_suffix, $eval, $np) = @_;

    # Create directory for sample
    my $assemlib = "$assemdir/$lib/";
    unless(-d $assemlib or mkdir $assemlib) { die "[ERROR assembleByProt $lib] Could not mkdir $assemlib\n"; }    

    # Set IO paths
    my $fil_1 = "$readdir/$lib" . "_$fwd_suffix.fastq.gz";
    my $fil_2 = "$readdir/$lib" . "_$rev_suffix.fastq.gz";
    my $fil_u = "$readdir/$lib" . "_$unpaired_suffix.fastq.gz";

    my $blast_1 = "$assemlib/$lib" . "_$fwd_suffix.against_targets.blast";
    my $blast_2 = "$assemlib/$lib" . "_$rev_suffix.against_targets.blast";
    my $blast_u = "$assemlib/$lib" . "_$unpaired_suffix.against_targets.blast";

    my $fasta_1 = "$assemlib/$lib" . "_$fwd_suffix.fasta";
    my $fasta_2 = "$assemlib/$lib" . "_$rev_suffix.fasta";
    my $fasta_u = "$assemlib/$lib" . "_$unpaired_suffix.fasta";

    # Blast reads against target proteins
    unless(-e "$blast_1") { blastProts($fil_1, $blast_1, $fasta_1, $blast_dbs_dir, $eval, $np); }
    unless(-e "$blast_2") { blastProts($fil_2, $blast_2, $fasta_2, $blast_dbs_dir, $eval, $np); }
    unless(-e "$blast_u") { blastProts($fil_u, $blast_u, $fasta_u, $blast_dbs_dir, $eval, $np); }

    # For each target protein that was hit by reads from a file, collate those reads
    my $call_1  = getbest($assemlib, $fasta_1, $blast_1, "1" , $target_seqs_list);
    my $call_2  = getbest($assemlib, $fasta_2, $blast_2, "2" , $target_seqs_list);
    my $call_u  = getbest($assemlib, $fasta_u, $blast_u, "u" , $target_seqs_list);

    # Collate those reads with the same sequence ID from the paired file
    my $call_1p = getbest($assemlib, $fasta_2, $blast_1, "1p", $target_seqs_list);
    my $call_2p = getbest($assemlib, $fasta_1, $blast_2, "2p", $target_seqs_list);
}

sub blastProts {
    my ($fil, $blast, $fasta, $blast_dbs_dir, $eval, $np) = @_;
    
    my $blast_call = "$blastall_path -p blastx -d $blast_dbs_dir/target_proteins -e $eval -m 8 -I T";
    my $start_time = localtime();
    print "Blasting $fil against target proteins with legacy BLAST at $start_time\n";

    # Read sample reads file, convert to FASTA and save a copy, split into blocks and BLASTx
    system(qq(
        zcat $fil | 
        awk '{if(NR % 4 == 1 || NR % 4 == 2) {sub(/@/, ">"); print; } }' | 
        tee $fasta |
        parallel -j $np --block 200K --recstart '>' --pipe $blast_call > $blast
    ));
}

sub getbest {
    my ($assemlib, $fasta, $blast, $f, $target_seqs_list) = @_;

    open TARGET_SEQS_LIST, "<$target_seqs_list" or die "[ERROR assembleByProt $lib] Could not open the target IDs list";
    my @protnames = <TARGET_SEQS_LIST>;
    chomp(@protnames);
    close(TARGET_SEQS_LIST);

    # Hash of sequence IDs -> sequences
    my %fasta;
    open(FA, "<$fasta") or die "[ERROR assembleByProt $lib] Failed to open FASTA file $fasta\n";
    while(<FA>) {
        # This is fine as read files have alternating ID read lines
        my $seqnm = $_;
        my $seq   = <FA>;
        $seqnm =~ s/>//;
        chomp($seqnm); chomp($seq);
        $fasta{ $seqnm } = $seq;
    }
    close FA;
        
    # Hash of protein names -> sequence IDs aligned to the prot -> sequences
    my %prothits = map { $_ => {} } @protnames;
    open(BLOUT, "<$blast") or die "[ERROR assembleByProt $lib] Failed to open BLAST output file $blast\n";
    while(<BLOUT>) {
       my @linbits = split(/\t/, $_);
       my $hitseqnm = $linbits[0];
       my $hitsprot = $linbits[1];
       $prothits{$hitsprot}{$hitseqnm} = $fasta{ $hitseqnm };
    }

    # For each protein that received hits
    foreach my $prot (keys %prothits) {
        # Create script output directory structure for target
        my $poutdir = "$assemlib/$prot/";
        unless(-d $poutdir or mkdir $poutdir) { die "[ERROR assembleByProt $lib] Could not create target protein output directory $poutdir\n"; }

        my $assemble_by_prot_dir = "$poutdir/${prot}_assemble_by_prot/";
        unless(-d $assemble_by_prot_dir or mkdir $assemble_by_prot_dir) {
            die "[ERROR assembleByProt $lib] Could not create assembleByProt output directory $assemble_by_prot_dir\n";
        }

        # Append the reads to the output file
        my $poutfil = "$assemble_by_prot_dir/${prot}_${f}_hitreads.fasta";
        open(POUT, ">$poutfil") or die "[ERROR assembleByProt $lib] Failed to open output hitreads file $poutfil\n"; 
        foreach my $seqhitprot (keys %{$prothits{$prot}}) {
            print POUT ">$seqhitprot\n$prothits{$prot}{$seqhitprot}\n";
        }
        close(POUT);
    }
}

__END__

=head1 NAME

assembleByProtv2 - BLASTx sample reads onto target proteins and collate reads that align from each file.

=head1 USAGE

=over

=item B<perl assembleByProtv2.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $readdir

Directory containing cleaned read files.

=item $assemdir

Output directory.

=item $fwd_suffix

Filename suffix for forward read files.

=item $rev_suffix

Filename suffix for reverse read files.

=item $unpaired_suffix

Filename suffix for unpaired read files.

=item $target_seqs_list

Text file with target protein IDs.

=item $np

Number of BLASTx processes to use.

=item $blastall_path

Path to blastall binary.

=item $eval

BLASTx expectation value.

=back

=head1 SUBROUTINES

=over

=item filtAssemb

Locate input filenames based on $*_suffix arguments and call blastProts and getbest.

=item blastProts

BLASTx sample reads $fil against target proteins database with expectation value of $eval using $np processes.

=item getbest

Collate reads from sample reads FASTA file $fasta that aligned in the BLASTx output file $blast to each target in $target_seqs_list.

=back

=head1 DIAGNOSTICS

=over 

=item [ERROR assembleByProt $lib] ... 

File or directory creation failed. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

NCBI blastall

GNU parallel

=head1 KNOWN BUGS

None.

=head1 NOTES

Uses legacy blastall due to performance gains over BLAST+ with short reads as queries.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
