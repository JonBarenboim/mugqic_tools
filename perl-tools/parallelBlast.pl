#!/usr/bin/env perl

=head1 NAME

I<parallelBlast>

=head1 SYNOPSIS

parallelBlast [options] [command [arguments]]

=head1 DESCRIPTION

B<parallelBlast> is a script made to do multiple BLAST
system calls in a multicore machine.

B<parallelBlast> will try to guess the number of cores 
on a linux machine, split the <query file> in 'number_of_core - 1' 
(using the program fastasplit) and  run in parallel one BLAST job 
per splitted <query file>.


B<parallelBlast> also accepts an specific number of jobs.
 
=head1 OPTIONS

B<--file or -f >I<query file>

=over 10

Fullpath to the multifasta file that will be querried against a database.

=back 


B<--OUT or -O >I<Final output file>

=over 10

This option specifies the name of the final output. In other words, after BLAST computes the alignments of each splitted
file B<ParalleBlast> concatenates the files into a final file named according to this option.

=back

B<--nthreads or -n >I<number of threads>

=over 10

Number of BLAST threads running simultaneosly. It will also determines in how many chunks the query file will be splitted.

=back

B<--BLAST or -b >I<BLAST options>

=over 10

I<IMPORTANT> the BLAST option must be between single quotes. 
q
Options to be passed to BLAST. Note that the BLAST option B<--query> should not be passed here. The script will pass this option to blast after
it splits the <query file>. The -o (output) option should not be passed direct to BLAST, since BLAST will output only partial results of 
the initial <query file>. Instead use the option B<--OUT or -O> explained above. I<SEE> example below.


=back



B<I<EXAMPLE>> Run BLASTx with 20 threads.

B<parallelBlast -file query_file.txt --OUT final_blast_results.txt -n 20 --BLAST 'blastx -db nr -outfmt 6 -evalue e-10'> 


=head1 AUTHOR

B<David Morais> - I<dmorais@cs.bris.ac.uk>

=head1 DEPENDENCY

B<Getopt::Long> Used to parse command line options.

B<Proc::ParallelLoop> parallel jobs 

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<File::Path> File path parsing

B<File::Basename> Path parsing

B<use Cwd> Path parsing

B<fastasplit> Split query_file
(www.bcgsc.ca/downloads/parts/software/resources/src/exonerate-1.4.0/src/util/fastasplit.c). 
=cut


# Strict Pragmas
#-------------------------- 

use strict;
use warnings;
#--------------------------

# INCLUDES
#--------------------------
use Pod::Usage;
use Data::Dumper;
use Proc::ParallelLoop;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd 'abs_path';
#--------------------------


# GLOBAL VARIABLES
#--------------------------
my $help;
my $BLAST;
my $outFile;
my $queryFile;
my $nthreads;
my $queryDir;
my $file;
my $splitfasta = dirname(abs_path($0));
#--------------------------



GetOptions(
  'help|h!'      => \$help,
  'file|f=s'     => \$queryFile,
  'OUT|O=s'      => \$outFile,
  'BLAST|b=s'    => \$BLAST,
  'nthreads|n=i' => \$nthreads
) or die "Fatal Error: Problem parsing command-line " . $!;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $help;



# Get Number of CPUs if nthreads was not specified
$nthreads = $nthreads ? $nthreads : `cat /proc/cpuinfo | grep "processor" | awk '{print \$NF}' | sort -u | wc -l`;
chomp $nthreads; 

#  MAIN
#--------------------

# Create temp dir
#--------------------
$queryDir = dirname($queryFile);
$file = basename($queryFile);
mkpath ("$queryDir/temp_$file") unless (-d "$queryDir/temp_$file");


# Split files
#---------------------
system "module load mugqic/exonerate && fastasplit -c $nthreads -f $queryFile -o $queryDir/temp_$file";

opendir(DIR, "$queryDir/temp_$file");
my @job = grep {/chunk/} (readdir(DIR));

# Start parallel jobs
#---------------------
pareach([@job], sub {
  my $sample = shift;
  
  # Parse BLAST command line
  #--------------------------
  $BLAST =~ s/-db/-query $queryDir\/temp\_$file\/$sample -db /;
  $BLAST .= " -out $queryDir/temp_$file/$sample.OUT";
  
  print $BLAST . "\n";
  system $BLAST;
  
 }, {"Max_Workers" => $nthreads});


# Concatenate files
#-------------------

opendir(DIR, "$queryDir/temp_$file");
my @blastResult = grep {/OUT/} (readdir(DIR));

foreach (@blastResult){
  system "cat $queryDir/temp_$file/$_ >>$outFile";
}

# Remove temp dir
#--------------------
system "rm -r $queryDir/temp_$file";
