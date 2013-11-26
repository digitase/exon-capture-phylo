use warnings;
use strict;

my $assemdir = "/home2/jgb/camaxlibs/";

# from the sequenced genome
# one anolis protein per fasta: more than 3k
# a list of the 3k proteins
my $protfil  = "/home2/jgb/assemble/target/anolisproteinlist.txt";

# protein sequences of targets
my $protseqs = "/home2/jgb/assemble/target/refs/anolistargetproteins/";


my $exonerate = "/usr/local/bin/exonerate";
my $cap3      = "/home/jgb/bin/cap3";

my @k            = (31, 41, 51, 61, 71, 81);

open PROTS, "<$protfil" or die "could not open the lib file";

my @prots = <PROTS>;
chomp(@prots); 
close(PROTS); 

my $lib = $ARGV[0];

foreach my $prot (@prots) {
   
   sleep(2);

   my $protseq = $protseqs . $prot . ".fasta";
   
   my $concatcontigs     = $assemdir . $lib . "/" . $prot . "_velvetsixk.fa";
   my $contigscap3ed     = $assemdir . $lib . "/" . $prot . "_velvetsixk.fa.cap3out";
   my $contigsexonerated = $assemdir . $lib . "/" . $prot . "_velvetsixk.fa.cap3out.exonerate";


   my $catcall = "";
   for my $k (@k) {
      my $kfil = $assemdir . $lib . "/" . $prot . "_velvetk" . $k . "/contigs.fa";
      $catcall = $catcall . " " . $kfil;
   }

   system("cat $catcall > $concatcontigs");

   # 99% identity in a 20 bp overlap
   # perfect overlap for short overlaps
   system("cap3 $concatcontigs -o 20 -p 99");

   # concat consensus and singleton contigs
   system("cat $concatcontigs.cap.contigs $concatcontigs.cap.singlets > $contigscap3ed");


   my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
   system("$exonerate --model protein2genome -q $protseq -t $contigscap3ed --showvulgar no --showalignment no --ryo $ryo > $contigsexonerated");
   
}


