# arrange best contigs by exon
use warnings;
use strict;

# my $assemdir = "/home2/jgb/camaxlibs/";
# my $libfil   = $assemdir . "camaxalllibs.txt";
# my $protfil  = "/home2/jgb/assemble/target/anolisproteinlist.txt";
# my $exonlist = "/home2/jgb/assemble/target/refs/targets/targetexons.txt.all";
# my $ctgnums  = $assemdir . "countallcontigs.txt";

my ($assemdir, $libfil, $protfil, $exonlist) = @ARGV;

my $gathercontigs_dir = "$assemdir/gathercontigs/";
unless(-e $gathercontigs_dir or mkdir $gathercontigs_dir) { die "could not make $gathercontigs_dir \n"; }    

my $alldir = "$gathercontigs_dir/all_contigs/";
unless(-e $alldir or mkdir $alldir) { die "could not make $alldir \n"; }    

my $bestdir = "$gathercontigs_dir/best_contigs/";
unless(-e $bestdir or mkdir $bestdir) { die "could not make $bestdir \n"; }    

open LIBS, "<$libfil" or die "could not open the lib file";
open EXONS, "<$exonlist" or die "could not open the lib file";
my @libs  = <LIBS>;
my @exons = <EXONS>;
chomp(@libs); chomp(@exons);
close(LIBS); close(EXONS);

my $ctgnums = "$gathercontigs_dir/countallcontigs.txt";

my %contignum; 
foreach my $exonfile (@exons) {

    if ($exonfile =~ /(ENS\S+)_(exon\S+)/) { 

        my $prot = $1; 
        my $exon = $1 . "_" . $2;   

        my $exonallcontigs = $alldir . $exon . "_all_contigs.fasta";
        open ALLCON, ">$exonallcontigs" or die "cannot open a contig file";

        my $exonbestcontigs = $bestdir . $exon . "_best_contigs.fasta";
        open BESTCON, ">$exonbestcontigs" or die "cannot open a contig file";

        foreach my $lib (@libs) {
            my $assemlib = "$assemdir/$lib/";
            my $bestcontig_distrib_dir = "$assemlib/$prot/${prot}_bestcontig_distrib/";

            my $contigfilall = "$bestcontig_distrib_dir/${exon}_velvet_contigs.cap3ed.exonerated.filtered.fasta";      
            my $contigfilbest = "$bestcontig_distrib_dir/${exon}_velvet_contigs.cap3ed.exonerated.filtered.best_contig.fasta";      

            if (-e $contigfilall) {

                open CTGA, "<$contigfilall" or die "could not open the contig file";
                my @confilalllines = <CTGA>;

                my $alltmpstring = join('',@confilalllines);
                my @tmplines = split(/\n>/,$alltmpstring);

                my $tmpnum = ( scalar @tmplines );
                $contignum{$exon}{$lib} = $tmpnum;

                foreach my $tmplin (@tmplines) {
                    my ($name, $tmpseq) = split(/\n/,$tmplin, 2);
                    $name =~ s/>//;
                    $name = $lib."_".$name;
                    $tmpseq =~ s/\n//g;
                    print ALLCON ">$name\n$tmpseq\n";

                }

            } else {
                $contignum{$exon}{$lib} = 0;
            }

           if(-e $contigfilbest) {

                open CTGB, "<$contigfilbest" or die "could not open the contig file";
                my @confilbestlines = <CTGB>;

                my $besttmpstring = join('',@confilbestlines);
                my @tmplines = split(/\n>/,$besttmpstring);

                foreach my $tmplin (@tmplines) {
                    my ($name, $tmpseq) = split(/\n/,$tmplin, 2);
                    $name =~ s/>//;
                    $name = $lib;
                    $tmpseq =~ s/\n//g;
                    print BESTCON ">$name\n$tmpseq\n";

                }

            }


        }

        close(ALLCON);
        close(BESTCON);
    }
}

open CONNUMS, ">$ctgnums" or die "cannot open contig number outfile $ctgnums";

my $p = 0;
foreach my $exon (keys %contignum){

     if ($p==0) {
        print CONNUMS "sample:"; 
          foreach my $lib (keys %{$contignum{$exon}}){
             print CONNUMS "\t$lib";
          }
        print CONNUMS "\n";

      }
     $p++;

     print CONNUMS "$exon";
     foreach my $lib (keys %{$contignum{$exon}}){
        print CONNUMS "\t$contignum{$exon}{$lib}";
     }
     print CONNUMS "\n";
}
close(CONNUMS);


