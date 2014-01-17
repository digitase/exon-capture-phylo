
use strict;
use warnings;

my ($assemdir, $targets_blast_db_name, $all_blast_db_name, $all_prot_seqs, $target_seqs_list) = @ARGV;

# Create blast database directory
my $blast_dbs_dir = "$assemdir/blast_dbs/";
unless(-d $blast_dbs_dir or mkdir $blast_dbs_dir) {
    die "Could not create blast database output directory $blast_dbs_dir\n";
}

my $logfile = "$assemdir/blast_dbs/prepareData.log";

# Copy in all prots file, using the first field as the ID
system("awk '{print \$1}' $all_prot_seqs > $blast_dbs_dir/all_proteins.fasta");

# Create the all proteins blast db
unless(-e "$blast_dbs_dir/$all_blast_db_name.pin") {
    system("makeblastdb -dbtype prot -in $blast_dbs_dir/all_proteins.fasta -out $blast_dbs_dir/$all_blast_db_name > $logfile");
}

# Extract sequences in target list from proteins file
my $target_seqs = "$blast_dbs_dir/target_proteins.fasta";
system("perl -ne 'if (/^>(\\S+)/) {\$c=\$i{\$1}}\$c?print:chomp;\$i{\$_}=1 if \@ARGV' $target_seqs_list $blast_dbs_dir/all_proteins.fasta > $target_seqs");

# Create the target BLAST database unless it already exists
unless(-e "$blast_dbs_dir/$targets_blast_db_name.pin") {
    system("makeblastdb -dbtype prot -in $target_seqs -out $blast_dbs_dir/$targets_blast_db_name > $logfile");
}

