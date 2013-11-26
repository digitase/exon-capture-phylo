use warnings;
use strict;

my $assemdir = "/home2/jgb/camaxlibs/";
my $exondir  = "/home2/jgb/assemble/crypto/refs/targets/";
my $exonlist = "/home2/jgb/assemble/crypto/refs/targets/targetexons.txt.all";
my $libfil   = $assemdir . "camaxalllibs.txt";
my $protfil  = "/home2/jgb/assemble/crypto/anolisproteinlist.txt";
my $protseqs = "/home2/jgb/assemble/crypto/refs/anolistargetproteins/";
my $blastdb  = "/home2/jgb/blastdb/anolis_carolinensis/Anolis_carolinensis.AnoCar2.0.67.pep.all.fa";

my $minoverlap = 0.65;

 open LIBS, "<$libfil" or die "could not open the lib file";
 open EXONS, "<$exonlist" or die "could not open the lib file";


 my @libs  = <LIBS>;
 my @exons = <EXONS>;
 chomp(@libs); chomp(@exons);
 close(LIBS); close(EXONS);

 foreach my $exonfile (@exons) {
sleep(2);
   if ($exonfile =~ /(ENS\S+)_(exon\d+)_/) { 

    my $prot = $1; 
    my $exon = $1 . "_" . $2;   

    my ($lower, $upper) = getlimits( $exondir.$exonfile );


#    foreach my $lib (@libs) {

      my $lib = $ARGV[0];
      my $contigsallkexonerate = $assemdir . $lib . "/" . $prot . "_velvetsixk.fa.cap3out.exonerate";
      my $exonlibfil = $assemdir . $lib . "/" . $exon . ".fa";      


      my $call1 = parseexon($contigsallkexonerate, $exonlibfil, $lower, $upper);
      my $call2 = performRBH($exonlibfil, $prot, $blastdb);



   }
   
   else {
      print "could not recognise the exon...";
   }



} # end foreach @exons


   

sub performRBH {
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
          if ($name eq $bestcontig) {
          $tmpseq =~ s/\n//g;
          print BEST ">$name\n$tmpseq\n";
       }

       }


    }

}


sub getlimits {

    my $exfilex = $_[0] . ".exonerate";
    open EFE, "<$exfilex" or die "could not open exonfile";
    my @lines = <EFE>;
    my $nameline = $lines[2];
    if ($nameline =~ / b(\d+) e(\d+) p/) {my $b = $1; my $e = $2; return($b, $e);}
    else {my $b = "0"; my $e = "0"; return($b, $e);}

}


sub parseexon {

    my $exonerated = $_[0];
    my $unexonerated = $_[1];
    my $lower = $_[2];
    my $upper = $_[3];

    open UNEX, ">$unexonerated" or die "could not open exonfile";
    open EX,   "<$exonerated" or die "could not open exonfile";

    #remove exonerate lines
    my @filtlines = ();
    my @lines = <EX>;
    foreach my $line (@lines) {
       if ($line =~ /^Command|^Host|completed exonerate|^\n$/) {
           # do nothing
       }
       else {
            push(@filtlines, $line);
       }

    }
    
    # check overlap with exon and print
    my $newfilstring = join('',@filtlines);
    my @tmplines = split(/\n>/,$newfilstring);

    my $t = 1;
    foreach my $tmplin (@tmplines) {
        my ($name, $tmpseq) = split(/\n/,$tmplin, 2);
        $name =~ s/>//;
        my $b; my $e; 
        if ($name =~ / b(\d+) e(\d+) p/) {$b = $1; $e = $2;}
        else {$b = "0"; $e = "0";}

        # is there any overlap 
        if ( ($b >= $lower && $b < $upper) || ($e > $lower && $e <= $upper) || ($b <= $lower && $e >= $upper)) {
  
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


  
