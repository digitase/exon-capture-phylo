use warnings;
use strict;

my ($lib, $assemdir, $k) = @ARGV;

use feature 'say';
say "Command line arguments to $0 are:";
say for @ARGV;

my $assemlib = $assemdir . $lib . "/";
opendir(DIR, $assemlib);

my @filprotnames;
my @files = readdir(DIR);

closedir(DIR);

foreach(@files) {
    if ($_ =~ /(\S+)_(u|1|2|1p|2p)_hitreads.fa/) {
        push(@filprotnames, $1);
    }
}

my %unique = map { $_ => 1 } @filprotnames;
my @protnames = keys %unique;

foreach my $protname (@protnames) {

    # slow down i/o
    #sleep(3); 

    my $assem_filu = $assemlib . $protname . "_u_hitreads.fa";
    my $assem_fil1 = $assemlib . $protname . "_1_hitreads.fa";
    my $assem_fil2 = $assemlib . $protname . "_2_hitreads.fa";
    my $assem_fil1p = $assemlib . $protname . "_1p_hitreads.fa";
    my $assem_fil2p = $assemlib . $protname . "_2p_hitreads.fa";


    my $filcall = "";
    if(-e $assem_filu) {$filcall = "$filcall $assem_filu";}
    if(-e $assem_fil1) {$filcall = "$filcall $assem_fil1";}
    if(-e $assem_fil2) {$filcall = "$filcall $assem_fil2";}
    if(-e $assem_fil1p) {$filcall = "$filcall $assem_fil1p";}
    if(-e $assem_fil2p) {$filcall = "$filcall $assem_fil2p";}

    my $veldir = $assemlib."/". $protname . "_velvetk$k/";
    unless(-e $veldir or mkdir $veldir)
    {die "could not make $veldir \n";}
    my $allhitreads  = $veldir . $protname . "_all_hitreads.fa";

    system("cat $filcall > $allhitreads");
    system("velveth $veldir $k -short -fasta $allhitreads");
    # some of these params are only appropriate for small assemblies
    system("velvetg $veldir -very_clean yes -max_branch_length 320 -max_gap_count 6 -cov_cutoff 5");

}
