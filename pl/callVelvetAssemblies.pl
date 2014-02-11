use warnings;
use strict;

my ($lib, $assemdir, $target_seqs_list, $velveth_path, $velvetg_path, $k) = @ARGV;

my $time = localtime();
print "callVelvetAssemblies on $lib with k = $k at $time\n";

my $assemlib = "$assemdir/$lib/";

open TARGET_SEQS_LIST, "<$target_seqs_list" or die "could not open the target IDs list";
my @protnames = <TARGET_SEQS_LIST>;
chomp(@protnames);
close(TARGET_SEQS_LIST);

foreach my $protname (@protnames) {
    my $veldir = "$assemlib/$protname/${protname}_call_velvet_assemblies/";
    unless(-e $veldir or mkdir $veldir) { die "Could not make $veldir\n"; }

    my $kvalue_dir = "$veldir/${protname}_k$k/";
    unless(-e $kvalue_dir or mkdir $kvalue_dir) { die "Could not make $kvalue_dir\n"; }

    my $allhitreads  = "$kvalue_dir/${protname}_all_hitreads.fasta";
    system("cat $assemlib/$protname/${protname}_assemble_by_prot/* > $allhitreads");

    # Paramters optimised for small assemblies
    system("$velveth_path $kvalue_dir $k -short -fasta $allhitreads > $kvalue_dir/${protname}_velveth.log");
    system("$velvetg_path $kvalue_dir -very_clean yes -max_branch_length 320 -max_gap_count 6 -cov_cutoff 5 > $kvalue_dir/${protname}_velvetg.log");

    my $all_assembled_contigs = "$veldir/${protname}_velvet_contigs.fasta";
    system("cat $kvalue_dir/contigs.fa >> $all_assembled_contigs");

    # sleep(3); 
}
