use warnings;
use strict;

my ($lib, $assemdir, $protfil, $protseqs, @k_values) = @ARGV;

open PROTS, "<$protfil" or die "could not open the lib file";
my @prots = <PROTS>;
chomp(@prots); 
close(PROTS);

foreach my $prot (@prots) {
   
   #sleep(2);
   my $protseq = $protseqs . $prot . ".fasta";
   
   # velvet contigs from all k values combined with cat
   my $concatcontigs     = $assemdir . $lib . "/$prot/" . $prot . "_velvetsixk.fa";

   # cap3 .contigs and .singlets files combined with cat
   my $contigscap3ed     = $assemdir . $lib . "/$prot/" . $prot . "_velvetsixk.fa.cap3out";

   # cap3ed contigs through exonerate
   my $contigsexonerated = $assemdir . $lib . "/$prot/" . $prot . "_velvetsixk.fa.cap3out.exonerate";

# cat velvet contigs for all k values into one file
   my $catcall = "";
   for my $k (@k_values) {
      my $kfil = $assemdir . $lib . "/$prot/" . $prot . "_velvet_k" . $k . "/contigs.fa";
      $catcall = $catcall . " " . $kfil;
   }
   system("cat $catcall > $concatcontigs");

   # 99% identity in a 20 bp overlap
   # perfect overlap for short overlaps
   # system("cap3 $concatcontigs -o 20 -p 99");
   # Segfaults if no seqs found
   system("cap3 $concatcontigs -o 20 -p 99 >> cap3.log");

   # concat consensus and singleton contigs
   system("cat $concatcontigs.cap.contigs $concatcontigs.cap.singlets > $contigscap3ed");

   my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
   system("exonerate --model protein2genome -q $protseq -t $contigscap3ed --showvulgar no --showalignment no --ryo $ryo > $contigsexonerated");
   
}


