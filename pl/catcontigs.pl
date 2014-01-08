use warnings;
use strict;

my ($lib, $assemdir, $target_seqs_list, $target_seqs, @k_values) = @ARGV;

open PROTS, "<$target_seqs_list" or die "could not open the lib file";
my @prots = <PROTS>;
chomp(@prots); 
close(PROTS);

my $assemlib = "$assemdir/$lib/";

foreach my $prot (@prots) {

    my $catcontigs_dir = "$assemlib/$prot/${prot}_catcontigs/";
    unless(-e $catcontigs_dir or mkdir $catcontigs_dir) { die "Could not make $catcontigs_dir\n"; }

    # velvet contigs from all k values combined with cat
    my $all_assembled_contigs = "$assemlib/$prot/${prot}_call_velvet_assemblies/${prot}_velvet_contigs.fasta";
    system("cp $all_assembled_contigs $catcontigs_dir");
    $all_assembled_contigs = "$catcontigs_dir/${prot}_velvet_contigs.fasta";

    # 99% identity in a 20 bp overlap
    # perfect overlap for short overlaps
    system("cap3 $all_assembled_contigs -o 20 -p 99");

    # concat consensus and singleton contigs
    my $cap3ed_contigs = "$catcontigs_dir/${prot}_velvet_contigs.cap3ed.fasta";
    system("cat $all_assembled_contigs.cap.contigs $all_assembled_contigs.cap.singlets > $cap3ed_contigs");

    my $exonerate_query = "$catcontigs_dir/$prot.fasta";
    system("perl -ne 'if(/^>(\\S+)/) { \$c = grep {/^\$1\$/} qw($prot) } print if \$c' $target_seqs > $exonerate_query");

    my $ryo = '">%ti b%qab e%qae p%pi\\n%tas\\n"';
    my $exonerated_contigs = "$catcontigs_dir/${prot}_velvet_contigs.cap3ed.exonerated.fasta";
    system("exonerate --model protein2genome --query $exonerate_query --target $cap3ed_contigs --showvulgar no --showalignment no --ryo $ryo > $exonerated_contigs");

    #sleep(2");

    # while(<>) {
        
        # if(/^>(\S+)/) {
            # $c = grep {/^$1$/} "ID";
        # }
        # print if $c;
    
    # }

}


