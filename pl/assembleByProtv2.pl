# blastx reads against targets
# reads with significant hits to each exon per file

use strict;
use warnings;

use File::Basename;

# Sample name, samples directory, output directory, FASTA with all the targets to make the blast db, BLAST database name
my ($lib, $readdir, $assemdir, $target_seqs, $adb) = @ARGV;

# Create blast database directory
my $blast_dbs_dir = $assemdir . "/blast_dbs/";
unless(-d $blast_dbs_dir or mkdir $blast_dbs_dir) {
    die "Could not create blast database output directory $blast_dbs_dir\n";
}

# Create the target BLAST database unless it already exists
chdir("$blast_dbs_dir") or die "Cannot chdir to $blast_dbs_dir\n";
unless(-e "$adb.pin") {
    system("makeblastdb -dbtype prot -in $target_seqs -out $adb");
}
chdir("$assemdir") or die "Cannot chdir to $assemdir\n";

# Expectation value for blastx
my $eval = "1e-9";
# Number of parallel blast processes to use
my $np = "16";
# Use blastall blastx instead of blast+ blastx?
my $use_legacy_blast = 1;

# First iteration
my $assem_iter_1 = filtAssemb($readdir, $assemdir, $blast_dbs_dir, $lib, $adb, $eval, $np, $use_legacy_blast);
# Continue to iterate?

sub filtAssemb {
    my ($readdir, $assemdir, $blast_dbs_dir, $lib, $adb, $eval, $np, $use_legacy_blast) = @_;

    # Create directory for sample
    my $assemlib = $assemdir . $lib . "/";
    unless(-d $assemlib or mkdir $assemlib) { die "Could not mkdir $assemlib\n"; }    
    # chdir($assemlib);

    # Paths to read files
    # This section determines what the suffixes of the input read files need to be
    # TODO consider generic suffixes    

    my $fil_1 = "$readdir/$lib" . "_1_final.fastq.gz";
    my $fil_2 = "$readdir/$lib" . "_2_final.fastq.gz";
    my $fil_u = "$readdir/$lib" . "_u_final.fastq.gz";

    my $blast_1 = "$assemlib/$lib" . "_1_final.blast";
    my $blast_2 = "$assemlib/$lib" . "_2_final.blast";
    my $blast_u = "$assemlib/$lib" . "_u_final.blast";

    my $fasta_1 = "$assemlib/$lib" . "_1_final.fasta";
    my $fasta_2 = "$assemlib/$lib" . "_2_final.fasta";
    my $fasta_u = "$assemlib/$lib" . "_u_final.fasta";

    # do the blasting
    unless(-e "$blast_1") { blastProts($fil_1, $blast_1, $fasta_1, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }
    unless(-e "$blast_2") { blastProts($fil_2, $blast_2, $fasta_2, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }
    unless(-e "$blast_u") { blastProts($fil_u, $blast_u, $fasta_u, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast); }

    # for each exon that was hit by reads from a file, collate those reads
    my $call_1  = getbest($assemlib, $fasta_1, $blast_1, "1" );
    my $call_2  = getbest($assemlib, $fasta_2, $blast_2, "2" );
    my $call_u  = getbest($assemlib, $fasta_u, $blast_u, "u" );

    # for each exon that was hit by reads from a file, collate those reads with the same sequence ID from the paired file
    my $call_1p = getbest($assemlib, $fasta_2, $blast_1, "1p");
    my $call_2p = getbest($assemlib, $fasta_1, $blast_2, "2p");
}

sub blastProts {
    my ($fil, $blast, $fasta, $blast_dbs_dir, $adb, $eval, $np, $use_legacy_blast) = @_;
    
    my $blast_call;
    if($use_legacy_blast) {
        print "Generating $blast with legacy BLAST\n";
        $blast_call = "blastall -p blastx -d $blast_dbs_dir/$adb -e $eval -m 8 -I T"
    } else {
        print "Generating $blast BLAST+\n";
        $blast_call = "blastx -db $blast_dbs_dir/$adb -out $blast -evalue $eval -outfmt 6 -show_gis"
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
    my ($assemlib, $fasta, $blast, $f) = @_;

    # Hash of sequence names -> sequences
    my %fasta;
    open(FA, "<$fasta") or die "Failed to open FASTA file $fasta\n";
    while (<FA>)
    {
        my $seqnm = $_;
        my $seq   = <FA>;
        $seqnm =~ s/>//;
        chomp($seqnm); chomp($seq);
        $fasta{ $seqnm } = $seq;
    }
    close FA;
        
    # Hash of protein names -> hit sequence names to the prot -> sequences
    my %prothits;
    open(BLOUT, "<$blast") or die "Failed to open BLAST output file $blast\n";
    while(<BLOUT>) 
    {
       my @linbits = split(/\t/, $_);
       my $hitseqnm = $linbits[0];
       my $hitsprot = $linbits[1];
       $prothits{$hitsprot}{$hitseqnm} = $fasta{ $hitseqnm };
    }

    # Gather sequences hit to each protein target
    foreach my $prot (keys %prothits) {
        my $poutdir = "$assemlib/$prot/";
        unless(-d $poutdir or mkdir $poutdir) {
            die "Could not create exon output directory $poutdir\n";
        }
        my $poutfil = $poutdir . $prot . "_" . $f . "_hitreads.fasta";

        open(POUT, ">$poutfil") or die "Failed to open poutfil $poutfil\n"; 
        foreach my $seqhitprot (keys %{$prothits{$prot}}) {
            print POUT ">$seqhitprot\n$prothits{$prot}{$seqhitprot}\n";
            # print "$seqhitprot, $prot\n";
        }
        close(POUT);
    }
}


__END__


