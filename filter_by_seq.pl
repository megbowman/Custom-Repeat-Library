#!/usr/bin/perl -w
# This script will read one fasta file and write to a new fasta file just the identifer of those sequences bigger than a specific length. 
# Megan Bowman
# 15 November 2013

use strict;
use Getopt::Long();
use Bio::Seq;
use Bio::SeqIO; 

my $usage = "\nUsage: $0 --fastafile <fasta file name>  --output <output file> \n";

my ($fastafile, $in, $output);

Getopt::Long::GetOptions('fastafile=s' => \$fastafile,
			 'output=s' => \$output,);     

if (!-e $fastafile) {
    die "$fastafile does not exist!\n";
}
   
if (-e $output) {
    die "$output already exists!\n";
}

my ($id, $seqobj, $length, $out);

$in  = Bio::SeqIO->new(-file => "$fastafile", -format => 'Fasta');
$out = Bio::SeqIO->new (-format => 'fasta', -file => ">$output");


while ($seqobj = $in ->next_seq()) {
  $id = $seqobj->display_id();
  $length = $seqobj->length();
  if ($length < 50) {
    next;
  }
  else {
    $out->write_seq($seqobj);    
  }
}

exit;

