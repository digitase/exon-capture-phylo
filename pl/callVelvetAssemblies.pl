use warnings;
use strict;

my ($lib, $assemdir, $k) = @ARGV;

my $time = localtime();
print "callVelvetAssemblies with k = $k at $time\n";

my $assemlib = $assemdir . $lib . "/";
# Get list of proteins that were hit
my @protnames = split("\n", `ls -p $assemlib | grep "/" | sed s'/.\$//'`);

foreach my $protname (@protnames) {
    my $veldir = "$assemlib/$protname/" . $protname . "_velvet_k$k/";
    unless(-e $veldir or mkdir $veldir) { die "Could not make $veldir\n"; }

    my $allhitreads  = $veldir . $protname . "_all_hitreads.fasta";
    system("ls -d $assemlib/$protname/* | grep hitreads.fasta | xargs cat > $allhitreads");

    # Paramters optimised for small assemblies
    system("velveth $veldir $k -short -fasta $allhitreads >> velveth.log");
    system("velvetg $veldir -very_clean yes -max_branch_length 320 -max_gap_count 6 -cov_cutoff 5 >> velvetg.log");

    # Slow down i/o
    # sleep(3); 
}
