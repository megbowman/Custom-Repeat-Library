#!/usr/bin/perl -w
 
# CRL_final_blast_parse.pl
# This parses the blastp output of elements against uniref. Elements with significant hits to genes are removed, along with 50 bp upstream and downstream of the blast hit. Remaining sequence that is less than 50 bp is removed completely. Output will be those elements with no significant blast hits and remaining sequence from elements with blast hits. 

# Megan Bowman
# 24 April 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $usage = "\nUsage: $0 --MITE_blast <path of the MITE blast result File> --LTR_blast <path of allLTR blast results>  --ModelUn_blast <path of modeler unknown blast results> --Model_blast <path of modelerID blast results> --MITE_seq <path to MITE.lib> --LTR_seq <path to allLTR.lib> --ModelUn_seq <path to ModelerUnknown.lib> --Model_seq <path to ModelerID.lib> --output1 <name of no hit elements file> --output2 <name of cleaned elements file>\n";

my ($MITE, $LTR, $ModelUn, $Model, $output1, $output2, $MITE_seq, $LTR_seq, $ModelUn_seq, $Model_seq, $h_start1, $h_start2, $h_start3, $h_start4, $h_end1, $h_end2, $h_end3, $h_end4);

GetOptions('MITE_blast=s' => \$MITE,
	   'LTR_blast=s' => \$LTR,
	   'ModelUn_blast=s' => \$ModelUn,
	   'Model_blast=s' => \$Model,
	   'output1=s' => \$output1,
	   'MITE_seq=s' => \$MITE_seq,
	   'LTR_seq=s' => \$LTR_seq,
	   'ModelUn_seq=s' => \$ModelUn_seq,
	   'Model_seq=s' => \$Model_seq,
	   'output2=s' =>  \$output2);

if (!defined $MITE || !defined $LTR || !defined $ModelUn || !defined $Model || !defined $output1 || !defined $output2 || !defined $MITE_seq || !defined $LTR_seq || !defined $ModelUn_seq || !defined $Model_seq) {
  die $usage;
}

if (!-e $MITE || !-e $LTR || !-e $ModelUn || !-e $Model || !-e $MITE_seq || !-e $LTR_seq || !-e $ModelUn_seq || !-e $Model_seq) {
    die $usage;
}

if (-e $output1) {
  die "Output file already exists!\n";
}

if (-e $output2) {
  die "Output file already exists!\n";
}


my ($seqobj1, $seqobj2, $seqobj3, $seqobj4, %sequencehash1, %sequencehash2, %sequencehash3, %sequencehash4, $MITE_IN, $LTR_IN, $ModelUn_IN, $Model_IN); 

$Model_IN = Bio::SearchIO->new ( -format => 'blast', -file => "$Model"); 


open OUT1, ">MITE_Passed_Elements.txt"; 
open OUT2, ">MITE_Cleaned_Elements.txt";
open OUT3, ">MITE_Tossed_Junk.txt";

$MITE_seq = Bio::SeqIO->new ( -format => 'fasta', -file => "$MITE_seq"); 
$MITE_IN = Bio::SearchIO->new ( -format => 'blast', -file => "$MITE"); 

my %seq_hash1;

my $num_read1 = 0;
while ($seqobj1 = $MITE_seq->next_seq()) { 
  ++$num_read1;
  my $id1 = $seqobj1->display_id();
  my $seq1 = $seqobj1->seq();
  if (exists($seq_hash1{$id1})) {
    if ($seq1 = $seq_hash1{$id1}) {
      next;
    }
    if ($seq1 =~ /$seq_hash1{$id1}/) {
      next;
    }
  }
  $sequencehash1{$id1} = $seqobj1;
  $seq_hash1{$id1} = $seq1;
}


my (%passed_hash1, %cleaned_hash1);

while (my $result1 = $MITE_IN ->next_result()) { 
  my $query_name1 = $result1->query_name();
  my $query_length1 = $result1->query_length();
  my $num_hits1 = $result1->num_hits();
  if ($num_hits1 == 0) { 
    if (exists $passed_hash1{$query_name1}) {
      next;
    }
    else { 
      $passed_hash1{$query_name1} = $seq_hash1{$query_name1};
      next;
    }
  }
  else { 
    while (my $hit1 = $result1->next_hit()) {
      my $hit_sig1 = $hit1->significance();
      if ($hit_sig1 > 1e-010) {
	if (exists $passed_hash1{$query_name1}) {
	  next;
	}
	else { 
	  $passed_hash1{$query_name1} = $seq_hash1{$query_name1};
	  next;
	}
      }
      else { 
	my @start_array1 = ();
	my @end_array1 = ();
	my @sorted_start1 = ();
	my @sorted_end1 = ();
	my $hit_name1 = $hit1->name();
	while (my $hsp1 = $hit1->next_hsp) {
	  my $hsp_start1 = $hsp1 ->start('query');
	  my $hsp_end1 = $hsp1 ->end('query');
	  if ($hsp_start1 > $hsp_end1) {
	    my $hsp_start_int1;
	    my $hsp_end_int1;
	    $hsp_start1 = $hsp_start_int1;
	    $hsp_end1 = $hsp_end_int1;
	    $hsp_end1 = $hsp_start_int1;
	    $hsp_start1 = $hsp_end_int1;
	  }	  
	  push (@start_array1, $hsp_start1);
	  push (@end_array1, $hsp_end1);
	}
	@sorted_start1 = sort { $a <=> $b } @start_array1; 
	if ($sorted_start1[0] - 49 < 0) {
	  $h_start1 = $sorted_start1[0];
	}
	if ($sorted_start1[0] - 49 >= 0) {
	  $h_start1 = $sorted_start1[0] - 49;
	}
	  @sorted_end1 = sort { $b <=> $a } @end_array1;
	if ($sorted_end1[0] + 49 > $query_length1) {
	  $h_end1 = $sorted_start1[0];
	}
	if ($sorted_end1[0] + 49 <= $query_length1) { 
	  $h_end1 = $sorted_start1[0] + 49;
	}
	foreach my $key1 (keys %sequencehash1) {
	  if ($key1 eq $query_name1) {
	    if ($h_start1 == 1) { 
	      if ($h_end1 != $query_length1) {
		if (!exists $cleaned_hash1{$query_name1}) {
		  my $part_two1 = $sequencehash1{$key1}->subseq($h_end1, $query_length1);
		  my $new_seq1 = $part_two1;
		  if (length($new_seq1 lt 50)) {
		    print OUT3 ">$query_name1\n";
		    $cleaned_hash1{$query_name1} = 1;
		  }
		  else {  
		    print OUT2 ">$query_name1\n";
		    print OUT2 "$new_seq1\n";
		    $cleaned_hash1{$query_name1} = 1;		    
		  }
		}
		next;
	      }
	      else  {
		print OUT3 ">$query_name1\n";	
		$cleaned_hash1{$query_name1} = 1;		    
		next;
	      }
	    }
	    else { 
	      if ($h_end1 = $query_length1) {
		if (!exists $cleaned_hash1{$query_name1}) {
		  my $part_one1 = $sequencehash1{$key1}->subseq(1, $h_start1);
		  my $new_seq1 = $part_one1;
		  if (length($new_seq1 lt 50)) {
		    print OUT3 ">$query_name1\n";
		    $cleaned_hash1{$query_name1} = 1;		    	    
		  }
		  else { #(length($new_seq1 gt 50)) 
		    print OUT2 ">$query_name1\n";
		    print OUT2 "$new_seq1\n";
		    $cleaned_hash1{$query_name1} = 1;		    	    
		  }
		} 
		next;
	      }
	    }
	    if ($h_end1 != $query_length1) { 
	      if (!exists $cleaned_hash1{$query_name1}) {
		my $part_one1 = $sequencehash1{$key1}->subseq(1, $h_start1);
		my $part_two1 = $sequencehash1{$key1}->subseq($h_end1, $query_length1);
		my $new_seq1 = $part_one1 . $part_two1;
		if (length($new_seq1 lt 50)) {
		  print OUT3 ">$query_name1\n";		    
		  $cleaned_hash1{$query_name1} = 1;		    	    
		}
		else { 
		  print OUT2 ">$query_name1\n";
		  print OUT2 "$new_seq1\n";
		  $cleaned_hash1{$query_name1} = 1;		    	    
		}
	      }
	      next;
	    } 
	  }
	}
      }
    }
  }
}




foreach my $kee1 (keys %passed_hash1) {
  print OUT1 ">$kee1\n";
  print OUT1 "$passed_hash1{$kee1}\n";
}

open OUT4, ">LTR_Passed_Elements.txt";
open OUT5, ">LTR_Cleaned_Elements.txt";
open OUT6, ">LTR_Tossed_Junk.txt";


$LTR_seq = Bio::SeqIO->new ( -format => 'fasta', -file => "$LTR_seq"); 
$LTR_IN = Bio::SearchIO->new ( -format => 'blast', -file => "$LTR"); 

my %seq_hash2;

my $num_read2 = 0;
while ($seqobj2 = $LTR_seq->next_seq()) { 
  ++$num_read2;
  my $id2 = $seqobj2->display_id();
  my $seq2 = $seqobj2->seq();
  if (exists($seq_hash2{$id2})) {
    if ($seq2 = $seq_hash2{$id2}) {
      next;
    }
    if ($seq2 =~ /$seq_hash2{$id2}/) {
      next;
    }
  }
  $sequencehash2{$id2} = $seqobj2;
  $seq_hash2{$id2} = $seq2;
}


my (%passed_hash2, %cleaned_hash2);

while (my $result2 = $LTR_IN ->next_result()) { 
  my $query_name2 = $result2->query_name();
  my $query_length2 = $result2->query_length();
  my $num_hits2 = $result2->num_hits();
  if ($num_hits2 == 0) { 
    if (exists $passed_hash2{$query_name2}) {
      next;
    }
    else { 
      $passed_hash2{$query_name2} = $seq_hash2{$query_name2};
      next;
    }
  }
  if ($num_hits2 > 0) {
    while (my $hit2 = $result2->next_hit()) {
      my $hit_sig2 = $hit2->significance();
      if ($hit_sig2 > 1e-010) {
	if (exists $passed_hash2{$query_name2}) {
	  next;
	}
	else { 
	  $passed_hash2{$query_name2} = $seq_hash2{$query_name2};
	  next;
	}
      }
      else {
	my @start_array2 = ();
	my @end_array2 = ();
	my @sorted_start2 = ();
	my @sorted_end2 = ();
	my $hit_name2 = $hit2->name();
	while (my $hsp2 = $hit2->next_hsp) {
	  my $hsp_start2 = $hsp2 ->start('query');
	  my $hsp_end2 = $hsp2 ->end('query');
	  if ($hsp_start2 > $hsp_end2) {
	    my $hsp_start_int2;
	    my $hsp_end_int2;
	    $hsp_start2 = $hsp_start_int2;
	    $hsp_end2 = $hsp_end_int2;
	    $hsp_end2 = $hsp_start_int2;
	    $hsp_start2 = $hsp_end_int2;
	  }	  
	  push (@start_array2, $hsp_start2);
	  push (@end_array2, $hsp_end2);
	}
	@sorted_start2 = sort { $a <=> $b } @start_array2;
	if ($sorted_start2[0] - 49 < 0) {
	  $h_start2 = $sorted_start2[0];
	}
	if ($sorted_start2[0] - 49 >= 0) {
	  $h_start2 = $sorted_start2[0] - 49;
	}
	  @sorted_end2 = sort { $b <=> $a } @end_array2;
	if ($sorted_end2[0] + 49 > $query_length2) {
	  $h_end2 = $sorted_start2[0];
	}
	if ($sorted_end2[0] + 49 <= $query_length2) { 
	  $h_end2 = $sorted_start2[0] + 49;
	}
	foreach my $key2 (keys %sequencehash2) {
	  if ($key2 eq $query_name2) {
	    if ($h_start2 == 1) { 
	      if ($h_end2 != $query_length2) {
		if (!exists $cleaned_hash2{$query_name2}) {
		  my $part_two2 = $sequencehash2{$key2}->subseq($h_end2, $query_length2);
		  my $new_seq2 = $part_two2;
		  if (length($new_seq2 lt 50)) {
		    print OUT6 ">$query_name2\n";
		    $cleaned_hash2{$query_name2} = 1;
		  }
		  else {
		    print OUT5 ">$query_name2\n";
		    print OUT5 "$new_seq2\n";
		    $cleaned_hash2{$query_name2} = 1;		    
		  }
		}
		next;
	      }
	      else { 
		print OUT6 ">$query_name2\n";	
		$cleaned_hash2{$query_name2} = 1;		    
		next;
	      }
	    }
	    if ($h_start2 > 1) {
	      if ($h_end2 = $query_length2) {
		if (!exists $cleaned_hash2{$query_name2}) {
		  my $part_one2 = $sequencehash2{$key2}->subseq(1, $h_start2);
		  my $new_seq2 = $part_one2;
		  if (length($new_seq2 lt 50)) {
		    print OUT6 ">$query_name2\n";
		    $cleaned_hash2{$query_name2} = 1;		    	    
		  }
		  else { 
		    print OUT5 ">$query_name2\n";
		    print OUT5 "$new_seq2\n";
		    $cleaned_hash2{$query_name2} = 1;		    	    
		  }
		}
		next;
	      }
	      else { 
		if (!exists $cleaned_hash2{$query_name2}) {
		  my $part_one2 = $sequencehash2{$key2}->subseq(1, $h_start2);
		  my $part_two2 = $sequencehash2{$key2}->subseq($h_end2, $query_length2);
		  my $new_seq2 = $part_one2 . $part_two2;
		  if (length($new_seq2 lt 50)) {
		    print OUT6 ">$query_name2\n";		    
		    $cleaned_hash2{$query_name2} = 1;		    	    
		  }
		  else { 
		    print OUT5 ">$query_name2\n";
		    print OUT5 "$new_seq2\n";
		    $cleaned_hash2{$query_name2} = 1;		    	    
		  }
		}
		next;
	      } 
	    }
	  }
	}
      }
    }
  }
}

foreach my $kee2 (keys %passed_hash2) {
  print OUT4 ">$kee2\n";
  print OUT4 "$passed_hash2{$kee2}\n";
}

 

$ModelUn_IN = Bio::SearchIO->new ( -format => 'blast', -file => "$ModelUn"); 
$ModelUn_seq = Bio::SeqIO->new ( -format => 'fasta', -file => "$ModelUn_seq"); 

open OUT7, ">ModelUn_Passed_Elements.txt"; 
open OUT8, ">ModelUn_Cleaned_Elements.txt";
open OUT9, ">ModelUn_Tossed_Junk.txt";

my %seq_hash3;

my $num_read3 = 0;
while ($seqobj3 = $ModelUn_seq->next_seq()) { 
  ++$num_read3;
  my $id3 = $seqobj3->display_id();
  my $seq3 = $seqobj3->seq();
  if (exists($seq_hash3{$id3})) {
    if ($seq3 = $seq_hash3{$id3}) {
      next;
    }
    if ($seq3 =~ /$seq_hash3{$id3}/) {
      next;
    }
  }
  $sequencehash3{$id3} = $seqobj3;
  $seq_hash3{$id3} = $seq3;
}


my (%passed_hash3, %cleaned_hash3);

while (my $result3 = $ModelUn_IN ->next_result()) { 
  my $query_name3 = $result3->query_name();
  my $query_length3 = $result3->query_length();
  my $num_hits3 = $result3->num_hits();
  if ($num_hits3 == 0) { 
    if (exists $passed_hash3{$query_name3}) {
      next;
    }
    else { 
      $passed_hash3{$query_name3} = $seq_hash3{$query_name3};
      next;
    }
  }
  if ($num_hits3 > 0) {
    while (my $hit3 = $result3->next_hit()) {
      my $hit_sig3 = $hit3->significance();
      if ($hit_sig3 > 1e-010) {
	if (exists $passed_hash3{$query_name3}) {
	  next;
	}
	else { 
	  $passed_hash3{$query_name3} = $seq_hash3{$query_name3};
	  next;
	}
      }
      else { 
	my @start_array3 = ();
	my @end_array3 = ();
	my @sorted_start3 = ();
	my @sorted_end3 = ();
	my $hit_name3 = $hit3->name();
	while (my $hsp3 = $hit3->next_hsp) {
	  my $hsp_start3 = $hsp3 ->start('query');
	  my $hsp_end3 = $hsp3 ->end('query');
	  if ($hsp_start3 > $hsp_end3) {
	    my $hsp_start_int3;
	    my $hsp_end_int3;
	    $hsp_start3 = $hsp_start_int3;
	    $hsp_end3 = $hsp_end_int3;
	    $hsp_end3 = $hsp_start_int3;
	    $hsp_start3 = $hsp_end_int3;
	  }	  
	  push (@start_array3, $hsp_start3);
	  push (@end_array3, $hsp_end3);
	}
	@sorted_start3 = sort { $a <=> $b } @start_array3;
	if ($sorted_start3[0] - 49 < 0) {
	  $h_start3 = $sorted_start3[0];
	}
	if ($sorted_start3[0] - 49 >= 0) {
	  $h_start3 = $sorted_start3[0] - 49;
	}
	  @sorted_end3 = sort { $b <=> $a } @end_array3;
	if ($sorted_end3[0] + 49 > $query_length3) {
	  $h_end3 = $sorted_start3[0];
	}
	if ($sorted_end3[0] + 49 <= $query_length3) { 
	  $h_end3 = $sorted_start3[0] + 49;
	}
	foreach my $key3 (keys %sequencehash3) {
	  if ($key3 eq $query_name3) {
	    if ($h_start3 == 1) { 
	      if ($h_end3 != $query_length3) {
		if (!exists $cleaned_hash3{$query_name3}) {
		  my $part_two3 = $sequencehash3{$key3}->subseq($h_end3, $query_length3);
		  my $new_seq3 = $part_two3;
		  if (length($new_seq3 lt 50)) {
		    print OUT9 ">$query_name3\n";
		    $cleaned_hash3{$query_name3} = 1;
		  }
		  else { 
		    print OUT8 ">$query_name3\n";
		    print OUT8 "$new_seq3\n";
		    $cleaned_hash3{$query_name3} = 1;		    
		  }
		}
		next;
	      }
	      else { 
		print OUT9 ">$query_name3\n";	
		$cleaned_hash3{$query_name3} = 1;		    
		next;
	      }
	    }
	    if ($h_start3 > 1) {
	      if ($h_end3 = $query_length3) {
		if (!exists $cleaned_hash3{$query_name3}) {
		  my $part_one3 = $sequencehash3{$key3}->subseq(1, $h_start3);
		  my $new_seq3 = $part_one3;
		  if (length($new_seq3 lt 50)) {
		    print OUT9 ">$query_name3\n";
		    $cleaned_hash3{$query_name3} = 1;		    	    
		  }
		  else { 
		    print OUT8 ">$query_name3\n";
		    print OUT8 "$new_seq3\n";
		    $cleaned_hash3{$query_name3} = 1;		    	    
		  }
		}
		next;
	      }
	      else { 
		if (!exists $cleaned_hash3{$query_name3}) {
		  my $part_one3 = $sequencehash3{$key3}->subseq(1, $h_start3);
		  my $part_two3 = $sequencehash3{$key3}->subseq($h_end3, $query_length3);
		  my $new_seq3 = $part_one3 . $part_two3;
		  if (length($new_seq3 lt 50)) {
		    print OUT9 ">$query_name3\n";		    
		    $cleaned_hash2{$query_name3} = 1;		    	    
		  }
		  else { 
		    print OUT8 ">$query_name3\n";
		    print OUT8 "$new_seq3\n";
		    $cleaned_hash3{$query_name3} = 1;		    	    
		  }
		}
		next;
	      } 
	    }
	  }
	}
      }
    }
  }
}



foreach my $kee3 (keys %passed_hash3) {
  print OUT7 ">$kee3\n";
  print OUT7 "$passed_hash3{$kee3}\n";
}


  
$Model_IN = Bio::SearchIO->new ( -format => 'blast', -file => "$Model"); 
$Model_seq = Bio::SeqIO->new ( -format => 'fasta', -file => "$Model_seq"); 

open OUT10, ">Model_Passed_Elements.txt"; 
open OUT11, ">Model_Cleaned_Elements.txt";
open OUT12, ">Model_Tossed_Junk.txt";

my %seq_hash4;

my $num_read4 = 0;
while ($seqobj4 = $Model_seq->next_seq()) { 
  ++$num_read4;
  my $id4 = $seqobj4->display_id();
  my $seq4 = $seqobj4->seq();
  if (exists($seq_hash4{$id4})) {
    if ($seq4 = $seq_hash4{$id4}) {
      next;
    }
    if ($seq4 =~ /$seq_hash4{$id4}/) {
      next;
    }
  }
  $sequencehash4{$id4} = $seqobj4;
  $seq_hash4{$id4} = $seq4;
}


my (%passed_hash4, %cleaned_hash4);

while (my $result4 = $Model_IN ->next_result()) { 
  my $query_name4 = $result4->query_name();
  my $query_length4 = $result4->query_length();
  my $num_hits4 = $result4->num_hits();
  if ($num_hits4 == 0) { 
    if (exists $passed_hash4{$query_name4}) {
      next;
    }
    else { 
      $passed_hash4{$query_name4} = $seq_hash4{$query_name4};
      next;
    }
  }
  if ($num_hits4 > 0) {
    while (my $hit4 = $result4->next_hit()) {
      my $hit_sig4 = $hit4->significance();
      if ($hit_sig4 > 1e-010) {
	if (exists $passed_hash4{$query_name4}) {
	  next;
	}
	else {
	  $passed_hash4{$query_name4} = $seq_hash4{$query_name4};
	  next;
	}
      }
      else { 
	my @start_array4 = ();
	my @end_array4 = ();
	my @sorted_start4 = ();
	my @sorted_end4 = ();
	my $hit_name4 = $hit4->name();
	while (my $hsp4 = $hit4->next_hsp) {
	  my $hsp_start4 = $hsp4 ->start('query');
	  my $hsp_end4 = $hsp4 ->end('query');
	  if ($hsp_start4 > $hsp_end4) {
	    my $hsp_start_int4;
	    my $hsp_end_int4;
	    $hsp_start4 = $hsp_start_int4;
	    $hsp_end4 = $hsp_end_int4;
	    $hsp_end4 = $hsp_start_int4;
	    $hsp_start4 = $hsp_end_int4;
	  }	  
	  push (@start_array4, $hsp_start4);
	  push (@end_array4, $hsp_end4);
	}
	@sorted_start4 = sort { $a <=> $b } @start_array4;
	if ($sorted_start4[0] - 49 < 0) {
	  $h_start4 = $sorted_start4[0];
	}
	if ($sorted_start4[0] - 49 >= 0) {
	  $h_start4 = $sorted_start4[0] - 49;
	}
	@sorted_end4 = sort { $b <=> $a } @end_array4;
	if ($sorted_end4[0] + 49 > $query_length4) {
	  $h_end4 = $sorted_start4[0];
	}
	if ($sorted_end4[0] + 49 <= $query_length4) { 
	  $h_end4 = $sorted_start4[0] + 49;
	}
	foreach my $key4 (keys %sequencehash4) {
	  if ($key4 eq $query_name4) {
 	    if ($h_start4 == 1) { 
	      if ($h_end4 != $query_length4) {
		if (!exists $cleaned_hash4{$query_name4}) {
		  my $part_two4 = $sequencehash4{$key4}->subseq($h_end4, $query_length4);
		  my $new_seq4 = $part_two4;
		  if (length($new_seq4 lt 50)) {
		    print OUT12 ">$query_name4\n";
		    $cleaned_hash4{$query_name4} = 1;
		  }
		  else {
		    print OUT11 ">$query_name4\n";
		    print OUT11 "$new_seq4\n";
		    $cleaned_hash4{$query_name4} = 1;		    
		  }
		}
		next;
	      }
	      if ($h_end4 = $query_length4) {
		print OUT12 ">$query_name4\n";	
		$cleaned_hash4{$query_name4} = 1;		    
		next;
	      }
	    }
	    if ($h_start4 > 1) {
	      if ($h_end4 = $query_length4) {
		if (!exists $cleaned_hash4{$query_name4}) {
		  my $part_one4 = $sequencehash4{$key4}->subseq(1, $h_start4);
		  my $new_seq4 = $part_one4;
		  if (length($new_seq4 lt 50)) {
		    print OUT12 ">$query_name4\n";
		    $cleaned_hash4{$query_name4} = 1;		    	    
		  }
		  else { #(length($new_seq4 gt 50)) { 
		    print OUT11 ">$query_name4\n";
		    print OUT11 "$new_seq4\n";
		    $cleaned_hash4{$query_name4} = 1;		    	    
		  }
		}
		next;
	      }
	      if ($h_end4 != $query_length4) {
		if (!exists $cleaned_hash4{$query_name4}) {
		  my $part_one4 = $sequencehash4{$key4}->subseq(1, $h_start4);
		  my $part_two4 = $sequencehash2{$key4}->subseq($h_end4, $query_length4);
		  my $new_seq4 = $part_one4 . $part_two4;
		  if (length($new_seq4 lt 50)) {
		    print OUT12 ">$query_name4\n";		    
		    $cleaned_hash4{$query_name4} = 1;		    	    
		  }
		  else { 
		    print OUT11 ">$query_name4\n";
		    print OUT11 "$new_seq4\n";
		    $cleaned_hash4{$query_name4} = 1;		    	    
		  }
		}
		next;
	      } 
	    }
	  }
	}
      }
    }
  }
}



foreach my $kee4 (keys %passed_hash4) {
  print OUT10 ">$kee4\n";
  print OUT10 "$passed_hash4{$kee4}\n";
}


system ("cat MITE_Passed_Elements.txt MITE_Cleaned_Elements.txt > final_MITE.lib");
system ("cat LTR_Passed_Elements.txt LTR_Cleaned_Elements.txt > final_LTR.lib");
system ("cat ModelUn_Passed_Elements.txt ModelUn_Cleaned_Elements.txt > final_ModelUn.lib");
system ("cat Model_Passed_Elements.txt Model_Cleaned_Elements.txt >  final_Model.lib");

system ("cat final_MITE.lib final_LTR.lib  final_Model.lib  final_ModelUn.lib > allRepeats.lib");
system ("cat final_MITE.lib final_LTR.lib final_Model.lib > KnownRepeats.lib");






























