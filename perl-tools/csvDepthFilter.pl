#!/usr/bin/env perl
#
use strict;
use String::Util 'trim';

my $inputFile = $ARGV[0];
my $minDepth = $ARGV[1];

open(FILE, $inputFile) or die "Can't open\n";
my $line = <FILE>;

my @values = split("\t",$line);
my @depthColumns;
for(my $i=0; $i < @values; $i++) {
  if($values[$i] =~ /Depth$/) {
    push(@depthColumns, $i);
  }
}
print STDOUT $line;
while($line = <FILE>) {
  my $doPrint = 0;
  @values = split("\t",$line);
  for my $i (@depthColumns) {
    my @depths = split(',', $values[$i]);
    my $no = trim($depths[0]);
    my $tu = trim($depths[1]);
    if($no > $minDepth || $tu > $minDepth) {
      $doPrint = 1;
      last;
    }
  }

  if($doPrint != 0) {
    print STDOUT $line;
  }
}
close(FILE);
