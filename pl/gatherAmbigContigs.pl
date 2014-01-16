# arrange best contigs by exon
use warnings;
use strict;

my ($assemdir, $libfil, $exonlist) = @ARGV;

my $gatherambigcontigs_dir = "$assemdir/gatherAmbigcontigs/";
unless(-e $gatherambigcontigs_dir or mkdir $gatherambigcontigs_dir) { die "could not make $gatherambigcontigs_dir \n"; }    

open LIBS, "<$libfil" or die "could not open the lib file";
my @libs  = <LIBS>;
close(LIBS);
chomp(@libs);

open EXONS, "<$exonlist" or die "could not open the exon file";
my @exons = <EXONS>;
close(EXONS);
chomp(@exons);

foreach my $exon_name (@exons) {
    
    my $exon_best_ambig_contigs = "$gatherambigcontigs_dir/${exon_name}_best_ambig_contig.fasta";

    foreach my $lib (@libs) {

        my $vcf2ambigfasta_refs_fil = "$assemdir/$lib/${lib}_vcf2ambigfasta/${lib}_best2refs.vcf2ambigfasta_refs.fasta";
        system(qq(
            perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($exon_name) } print if \$c' $vcf2ambigfasta_refs_fil |
            awk '/^>/ {print \$0} !/^>/ {printf "%s", \$0} END {print}' >> $exon_best_ambig_contigs
        ));

    }

}

