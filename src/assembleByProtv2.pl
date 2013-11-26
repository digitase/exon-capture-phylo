# blastx reads against targets
# reads with significant hits to each exon per file

use warnings;
use strict;


# nominate directories

  my $readdir  = "/home2/jgb/reads/";
  my $assemdir = "/home2/jgb/camaxlibs/";

# make the blast database 
# unless it already exists

  my $adb    = "/home2/jgb/assemble/target/refs/anolistargetproteins.fa";
  unless(-e "$adb.pin")
  {
     system("formatdb -i $adb -p T");
  }

# set a couple of parameters 
  my $eval = "1e-9";
  my $np = "1";

#my $libfil = $assemdir . "libraries.txt";
#open(LIB, "<$libfil");
#my @libs = <LIB>;
#chomp(@libs);

# main loop, over libraries
# foreach my $lib (@libs) {

    # library name
   my $lib = $ARGV[0];
   my $assem_iter_1 = filtAssemb($readdir, $assemdir, $lib, $adb, $eval, $np);
# }


sub filtAssemb {
   
   my ($readdir, $assemdir, $lib, $adb, $eval, $np) = @_;
   my $assemlib = $assemdir . $lib . "/";

   unless(-e $assemlib or mkdir $assemlib)
   {die "could not make $assemlib \n";}

   # grab files, set up   
   my $fil_1 = $lib . "_1_final.fastq";
   my $fil_2 = $lib . "_2_final.fastq";
   my $fil_u = $lib . "_u_final.fastq";

   system("cp $readdir$fil_1.gz $assemlib"); system("gunzip $assemlib$fil_1.gz");
   system("cp $readdir$fil_2.gz $assemlib"); system("gunzip $assemlib$fil_2.gz");
   system("cp $readdir$fil_u.gz $assemlib"); system("gunzip $assemlib$fil_u.gz");

   my $fastq_1 = $assemlib.$fil_1; my $fastq_2 = $assemlib.$fil_2; my $fastq_u = $assemlib.$fil_u;
   my $fasta_1 = $assemlib.$fil_1; my $fasta_2 = $assemlib.$fil_2; my $fasta_u = $assemlib.$fil_u;
   my $blast_1 = $assemlib.$fil_1; my $blast_2 = $assemlib.$fil_2; my $blast_u = $assemlib.$fil_u;

   $fasta_1 =~ s/fastq/fasta/; $fasta_2 =~ s/fastq/fasta/; $fasta_u =~ s/fastq/fasta/;
   $blast_1 =~ s/fastq/blast/; $blast_2 =~ s/fastq/blast/; $blast_u =~ s/fastq/blast/;

   # do the blasting, put hits in fastas, arranged by Anolis protein
   my $f;
   unless (-e $blast_1 ) { $f = "1"; my $call_1 = blastProts($assemlib, $fastq_1, $fasta_1, $blast_1, $adb, $eval, $np, $f); }
   unless (-e $blast_2 ) { $f = "2"; my $call_2 = blastProts($assemlib, $fastq_2, $fasta_2, $blast_2, $adb, $eval, $np, $f); }
   unless (-e $blast_u ) { $f = "u"; my $call_u = blastProts($assemlib, $fastq_u, $fasta_u, $blast_u, $adb, $eval, $np, $f); }

   $f = "u";  my $call_u  = getbest($assemlib, $fasta_u, $blast_u, $f);
   $f = "1";  my $call_1  = getbest($assemlib, $fasta_1, $blast_1, $f);
   $f = "2";  my $call_2  = getbest($assemlib, $fasta_2, $blast_2, $f);
   $f = "1p"; my $call_1p = getbest($assemlib, $fasta_2, $blast_1, $f);
   $f = "2p"; my $call_2p = getbest($assemlib, $fasta_1, $blast_2, $f);

}


sub blastProts{

    my ($assemlib, $fastq, $fasta, $blast, $adb, $eval, $np, $f) = @_;

    open(FQ, "<$fastq");
    open(FA, ">$fasta");

    my %fasta;
    while (<FQ>)
    {
        my $seqnm = $_;
        my $seq   = <FQ>;
        my $qnam  = <FQ>;
        my $qscor = <FQ>;
        $seqnm =~ s/@//;
        chomp($seqnm); chomp($seq);
        print FA ">" . $seqnm ."\n". $seq . "\n";
        $fasta{ $seqnm } = $seq;
        #print "$fasta{ $seqnm }";
    }
    close FQ;
    close FA;

    my $callblast = system("blastall -p blastx -i $fasta -d $adb -o $blast -e $eval -m 8 -a $np -I T");
}



sub getbest{

    my ($assemlib, $fasta, $blast, $f) = @_;

    open(FA, "<$fasta");

    my %fasta;
    while (<FA>)
    {
        my $seqnm = $_;
        my $seq   = <FA>;
        $seqnm =~ s/>//;
        chomp($seqnm); chomp($seq);
        $fasta{ $seqnm } = $seq;
     }
     close FA;
        

    my %prothits;

    open(BLOUT, "<$blast");
    while (<BLOUT>)
    {
       my @linbits = split(/\t/, $_);
       my $hitseqnm = $linbits[0];
       my $hitsprot = $linbits[1];
       #print "$hitseqnm \t $hitsprot\n";
       $prothits{$hitsprot}{$hitseqnm} = $fasta{ $hitseqnm };
       #print "$fasta{ $hitseqnm }\n"
    }

# loop over protein targets
    foreach my $prot (keys %prothits){
    
        my $poutfil = $assemlib . $prot . "_" . $f . "_hitreads.fa";
        open(POUT, ">$poutfil"); 
          foreach my $seqhitprot (keys %{$prothits{$prot}}){
            print POUT ">$seqhitprot\n$prothits{$prot}{$seqhitprot}\n";
          }
        close(POUT);
      
    }
}

__END__


