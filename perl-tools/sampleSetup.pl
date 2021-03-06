#!/usr/bin/env perl

use Cwd;
use File::Basename;
use File::Path qw(mkpath);
use Getopt::Long;
#use Text::CSV;
use Text::CSV::Encoded;
use File::Find;

my $version = "1.0";

use strict;

&main();

sub getUsage {
  my $usage = <<END;
Usage: perl $0 --nanuqAuthFile \$HOME/.nanuqAuth.txt --usesheet project.nanuq.csv --tech HiSeq
  --nanuqAuthFile  <FILE>          Path to Nanuq authentication file
  --projectId      <INT>           Nanuq project ID from which to get the sample sheet (can't be used with --usesheet)
  --usesheet       <FILE>          Use the specified sample sheet instead of fecthing it from Nanuq (can't be used with --projectId)
  --links                          Create raw_reads directory and symlinks (default)
  --nolinks                        Do not create raw_reads directory or symlinks
  --tech           [HiSeq|MiSeq]   Sequencing technology ('MiSeq' or 'HiSeq' or '454')
  --excludeSent                    Exclude samples that are flagged as previously sent to client in nanuq
  --help                           Show this help

  The 'nanuqAuthFile' contains your Nanuq username and password.
  To create it, type the command:
  echo -n "user=<USERNAME>&password=<PASSWD>" > \$HOME/.nanuqAuth.txt ; chmod u+r,go-rwx \$HOME/.nanuqAuth.txt
  The '-n' is important because there cannot be a Carriage Return, EOL at the end of the line.
END

  return $usage;
}

sub main {
  my $techName;
  my $projectId;
  my $sampleSheet;
  my $nanuqAuthFile;
  my $links = 1;  # Create symlinks by default
  my $excludeSent = 0;
  my $help;
  my $result = GetOptions(
    "tech=s"          => \$techName,
    "projectId=i"     => \$projectId,
    "usesheet=s"      => \$sampleSheet,
    "nanuqAuthFile=s" => \$nanuqAuthFile,
    "links!"          => \$links,
    "excludeSent!"    => \$excludeSent,
    "help!"           => \$help,
  );

  if ($help) {
    die getUsage();
  }

  my $errMsg = "";
  if (defined($projectId) and defined($sampleSheet)) {
    $errMsg .= "Error: --projectId and --useSheet options cannot be used together!\n";
  }
  if (defined($projectId) and (!defined($nanuqAuthFile) or !(-e $nanuqAuthFile))) {
    $errMsg .= "Error: missing nanuqAuthFile!\n";
  }
  if ((!defined($projectId) or length($projectId) == 0) and (!defined($sampleSheet) or length($sampleSheet) == 0)) {
    $errMsg .= "Error: missing --projectId or --useSheet option!\n";
  }
  if (!defined($techName) or length($techName) == 0 or ($techName ne "HiSeq" and $techName ne "MiSeq" and $techName ne "454")) {
    $errMsg .= "Error: missing or invalid --tech value (should be 'HiSeq' or 'MiSeq' or '454')!\n";
  }
  if (length($errMsg)) {
    die $errMsg . "\n" . getUsage();
  }

  # Default Nanuq project file name
  my $projectFile = 'project.nanuq.csv';

  if (defined($projectId)) {
    # Fecth sample sheet data from Nanuq
    getSampleSheet($projectFile, $techName, $projectId, $excludeSent, $nanuqAuthFile);
  } else {
    # Use the specified sample sheet
    $projectFile = $sampleSheet;
  }

  if ($links) {
    my $rA_sampleInfos = parseSampleSheet($projectFile, $techName);
    createLinks($rA_sampleInfos, $techName);
    downloadBEDs($rA_sampleInfos, $nanuqAuthFile);
  }
}

sub getSampleSheet {
  my $projectFile = shift;
  my $techName = shift;
  my $projectId = shift;
  my $excludeSent = shift;
  my $nanuqAuthFile = shift;

  my $command = 'wget --no-cookies --post-file ' . $nanuqAuthFile . ' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/' . $techName . '/project/' . $projectId;
  if($excludeSent) {
    $command .= '/excludeStatusList/resultSentToClient,analysing';
  }
  $command .= '/filename/' . $projectFile . "\n";

  print '#' . $command;
  system($command);
  if ($? == -1) {
    print "failed to execute: $!\n";
    exit(1);
  } elsif ($? & 127) {
    printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
    exit(1);
  } else {
    my $childValue = $? >> 8;
    if ($childValue != 0) {
      printf "child exited with value %d\n", $childValue;
      exit(1);
    }
  }
}

sub downloadBEDs {
  my $rA_sampleInfos = shift;
  my $nanuqAuthFile = shift;

  my %bedFiles;
  for my $rH_sample (@$rA_sampleInfos) {
    if(defined($rH_sample->{'bed'}) && length($rH_sample->{'bed'}) > 0) {
      my @bedFileList = split(';', $rH_sample->{'bed'});
      for my $bed (@bedFileList) {
        $bedFiles{$bed} = 1;
      }
    }
  }

  for my $bed (keys(%bedFiles)) {
    my $command = 'wget --no-cookies --post-file ' . $nanuqAuthFile . ' https://genomequebec.mcgill.ca/nanuqLimsCgi/targetRegion/downloadBed.cgi?bedName=' . $bed . ' -O ' . $bed;
    print '#' . $command;
    system($command);
    if ($? == -1) {
      print "failed to execute: $!\n";
      exit(1);
    } elsif ($? & 127) {
      printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
      exit(1);
    } else {
      my $childValue = $? >> 8;
      if ($childValue != 0) {
        printf "child exited with value %d\n", $childValue;
        exit(1);
      }
    }
  }
}

sub parseSampleSheet {
  my $fileName = shift;
  my $techName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $qualOffsetIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $readSetIdIdx=-1;
  my $filePrefixIdx=-1;
  my $fastq1Idx=-1;
  my $fastq2Idx=-1;
  my $bamIdx=-1;
  my $bedIdx=-1;
  my $dateIdx="";

  my $csv = Text::CSV::Encoded->new ({ encoding => "iso-8859-1" });
  $csv->parse($line);
  my @headers = $csv->fields();
  for (my $idx = 0; $idx < @headers; $idx++) {

    $headers[$idx] =~ s/"//g;
    if ($headers[$idx] eq "Name") {
      $nameIdx = $idx;
    } elsif ($headers[$idx] eq "Library Barcode") {
      $libraryBarcodeIdx = $idx;
    } elsif ($headers[$idx] eq "Run") {
      $runIdIdx = $idx;
    } elsif ($headers[$idx] eq "Quality Offset") {
      $qualOffsetIdx = $idx;
    } elsif ($headers[$idx] eq "Region") {
      $laneIdx = $idx;
    } elsif ($headers[$idx] eq "Run Type") {
      $runTypeIdx = $idx;
    } elsif ($headers[$idx] eq "Status") {
      $statusIdx = $idx;
    } elsif ($headers[$idx] eq "Read Set Id") {
      $readSetIdIdx = $idx;
    } elsif ($headers[$idx] eq "Filename Prefix") {
      $filePrefixIdx = $idx;
    } elsif ($headers[$idx] eq "FASTQ1") {
      $fastq1Idx = $idx;
    } elsif ($headers[$idx] eq "FASTQ2") {
      $fastq2Idx = $idx;
    } elsif ($headers[$idx] eq "BED Files") {
      $bedIdx = $idx;
    } elsif ($headers[$idx] eq "BAM") {
      $bamIdx = $idx;
    
    # 454:
    } elsif ($headers[$idx] eq "Run Start Date" && $techName eq "454") {
      $dateIdx = $idx;
    } #elsif ($headers[$idx] eq "Run Start Date") && $techName eq "454" {
    #  $dateIdx = $idx;
    #} elsif ($headers[$idx] eq "Run Start Date") && $techName eq "454" {
    #  $dateIdx = $idx;
    #}


  }

  my $sampleSheetErrors = "";
  my $sampleSheetWarnings = "";
  if ($nameIdx == -1) {
    $sampleSheetErrors .= "Missing Sample Name\n";
  }
  if ($libraryBarcodeIdx == -1) {
    $sampleSheetErrors .= "Missing Library Barcode\n";
  }
  if ($runIdIdx == -1) {
    $sampleSheetErrors .= "Missing Run ID\n";
  }
  if ($qualOffsetIdx == -1) {
    $sampleSheetErrors .= "Missing Quality Offset\n";
  }
  if ($laneIdx == -1) {
    $sampleSheetErrors .= "Missing Lane\n";
  }
  if ($runTypeIdx == -1) {
    $sampleSheetErrors .= "Missing Run Type\n";
  }
  if ($statusIdx == -1) {
    $sampleSheetErrors .= "Missing Status\n";
  }
  if ($readSetIdIdx == -1) {
    $sampleSheetErrors .= "Missing Read Set Id\n";
  }
  if ($filePrefixIdx == -1) {
    $sampleSheetErrors .= "Missing Filename Prefix\n";
  }
  if ($fastq1Idx == -1 and $bamIdx == -1) {
    $sampleSheetErrors .= "Missing FASTQ1 or BAM\n";
  }
  if ($bedIdx == -1) {
    $sampleSheetWarnings .= "Missing BED Files column\n";
  }
  if (length($sampleSheetWarnings) > 0) {
    warn $sampleSheetWarnings;
  }
  if (length($sampleSheetErrors) > 0 && $techName ne "454") {
    die $sampleSheetErrors;
  }

  while ($line = <SAMPLE_SHEET>) {
    $csv->parse($line);
    my @values = $csv->fields();
    if ($values[$statusIdx] ne 'Data is valid') {
      warn "[Warning] Sample Name $values[$nameIdx], Run ID $values[$runIdIdx], Lane $values[$laneIdx] data is not in valid state!\n";
    } else {
      my %sampleInfo;
      $sampleInfo{'name'}           = $values[$nameIdx];
      $sampleInfo{'libraryBarcode'} = $values[$libraryBarcodeIdx];
      $sampleInfo{'runId'}          = $values[$runIdIdx];
      $sampleInfo{'qualOffset'}     = $values[$qualOffsetIdx];
      $sampleInfo{'lane'}           = $values[$laneIdx];
      $sampleInfo{'runType'}        = $values[$runTypeIdx];
      $sampleInfo{'readSetId'}      = $values[$readSetIdIdx];
      $sampleInfo{'filePrefix'}     = $values[$filePrefixIdx];
      $sampleInfo{'date'}           = $values[$dateIdx];
      $sampleInfo{'bed'}            = $values[$bedIdx];

      my $rootDir;
      if ($techName eq 'HiSeq') {
        $rootDir = "/lb/robot/hiSeqSequencer/hiSeqRuns/";
      } elsif ($techName eq 'MiSeq') {
        $rootDir = "/lb/robot/miSeqSequencer/miSeqRuns/";
      } elsif ($techName eq '454') {
        $rootDir = "/lb/robot/454sequencers/runs/";
      }else {
        die "Unknown prefix technology type: " . $techName . "\n";
      }

      # If Illumina:
      if($techName ne "454"){
        if ($values[$fastq1Idx]) {
          $sampleInfo{'fastq1'} = $rootDir . $values[$fastq1Idx];
        }
        if ($values[$fastq2Idx]) {
          $sampleInfo{'fastq2'} = $rootDir . $values[$fastq2Idx];
        }
        if ($values[$bamIdx]) {
          $sampleInfo{'bam'} = $rootDir . $values[$bamIdx];
        }

        # Readsets must have at least one BAM or FASTQ file defined to be selected, otherwise warning is raised
        if ($values[$bamIdx] or $values[$fastq1Idx]) {
          push(@retVal, \%sampleInfo);
        } else {
          warn "[Warning] Sample Name $values[$nameIdx], Run ID $values[$runIdIdx], Lane $values[$laneIdx] has neither BAM nor FASTQ1 fields set!\n";
        }
      }else{
        # else if 454
        my $fastaName        = $values[$filePrefixIdx].".fna";
        $fastaName =~ s/\.fna/\.454Reads\.fna/;
        my $qualName         = $values[$filePrefixIdx].".qual";
        $qualName =~ s/\.qual/\.454Reads\.qual/;
        $sampleInfo{'fasta'} = $fastaName;
        $sampleInfo{'qual'}  = $qualName;
        push(@retVal, \%sampleInfo);
      }
    }
  }

  return \@retVal;
}

sub createLinks {
  my $rA_sampleInfos = shift;
  my $techName       = shift;

  my @symlinks;
  
  ###################
  ## 454 
 
  if($techName eq "454"){
    my $rawReadDir = 'raw_reads/';
    mkpath($rawReadDir);
    # Create base raw read directory

    my @searchDirs= ("/lb/robot/454sequencer/runs/", "/lb/robot/454sequencer/runs/2014/","/lb/robot/454sequencer/runs/2015/","/lb/robot/454sequencer/runs/2016/");
    my $rootDir;
    my @directories;
    for my $indir (@searchDirs) {
      opendir(D, $indir) || warn "Can't open directory: $!\n";
      while (my $f = readdir(D)) {  
        #print "\$f = $f\n";
        push(@directories, $indir.$f);
      }
      closedir(D);
    }
    for my $rH_sample (@$rA_sampleInfos) {

      # find directory. For 454 libraries, libraries are stored by their date name.
      my $completeDate = $rH_sample->{'date'};
      my @completeDate = split(/\s/, $completeDate);
      my $date         = $completeDate[0];
      my $time         = $completeDate[1];
      my @date         = split(/-/, $date);
      my @time         = split(/:/, $time);
      my $year         = $date[0];
      my $month        = $date[1];
      my $day          = $date[2];
      my $hour         = $time[0];
      my $min          = $time[1];
      my $sec          = $time[2];
      
      my $prefix = $year."_".$month."_".$day."_".$hour."_".$min."_".$sec."_FLX";

      my $runDir = undef;
      for my $directory (@directories) {
        if($directory =~ /$prefix/){
          $runDir = $directory.'/';
          last;
        }
      }
      if(!defined($runDir)) {
        die "Run dir not found for: ".$prefix."\n";
      }

      my $runDDir = undef;
      opendir(D, $runDir) || die "Can't open directory: $!\n";
      while (my $f = readdir(D)) {  
        #print "\$f = $f\n";
        if($f =~ "^D_"){
          $runDDir = $runDir.$f."/";
          last;
        }
      }
      closedir(D);
      if(!defined($runDDir)) {
        die "Run D dir not found in: ".$runDir."\n";
      }

      push(@symlinks, [$runDDir."/".$rH_sample->{'fasta'}, $rawReadDir."/".$rH_sample->{'filePrefix'}.".fna"]);
      push(@symlinks, [$runDDir."/".$rH_sample->{'qual'}, $rawReadDir."/".$rH_sample->{'filePrefix'}.".qual"]);
    } 

    # Create all symbolic links
    for my $symlink (@symlinks) {
      if (-l @$symlink[1]) {
        warn "[Warning] Symbolic link @$symlink[1] already exists! Skipping.\n";
      } elsif (-f @$symlink[0] and symlink(@$symlink[0], @$symlink[1])) {
        print "Created symbolic link @$symlink[1] successfully.\n";
      } else {
        die "[Error] Can't create symbolic link @$symlink[1] to target @$symlink[0]!\n";
      }
    }

  ###################
  # HiSeq and MiSeq
  }else{

    for my $rH_sample (@$rA_sampleInfos) {

      # Create base raw read directory
      my $rawReadDir = 'raw_reads/' . $rH_sample->{'name'} . "/run" . $rH_sample->{'runId'} . "_" . $rH_sample->{'lane'};
      mkpath($rawReadDir);

      my $rawReadPrefix = $rawReadDir . '/' . $rH_sample->{'name'} . '.' . $rH_sample->{'libraryBarcode'} . '.' . $rH_sample->{'qualOffset'} . ".";
  

     # List all links to create
     if ($rH_sample->{'bam'}) {
       push(@symlinks, [$rH_sample->{'bam'}, $rawReadPrefix . "bam"]);
      } elsif ($rH_sample->{'fastq1'}) {
       my $runType = $rH_sample->{'runType'};
       if ($runType eq "SINGLE_END") {
         push(@symlinks, [$rH_sample->{'fastq1'}, $rawReadPrefix . "single.fastq.gz"]);
       } elsif ($runType eq "PAIRED_END") {
         push(@symlinks, [$rH_sample->{'fastq1'}, $rawReadPrefix . "pair1.fastq.gz"]);
         push(@symlinks, [$rH_sample->{'fastq2'}, $rawReadPrefix . "pair2.fastq.gz"]);
       } else {
         die "Error: unknown run type: $runType!";
       }
     }
    }

    # Create all symbolic links
    for my $symlink (@symlinks) {
      if (-l @$symlink[1]) {
        warn "[Warning] Symbolic link @$symlink[1] already exists! Skipping.\n";
      } elsif (-f @$symlink[0] and symlink(@$symlink[0], @$symlink[1])) {
        print "Created symbolic link @$symlink[1] successfully.\n";
      } else {
        die "[Error] Can't create symbolic link @$symlink[1] to target @$symlink[0]!\n";
      }
    }
  }
}
