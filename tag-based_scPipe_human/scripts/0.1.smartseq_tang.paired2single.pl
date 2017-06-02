#!/usr/bin/env perl
use strict;
use warnings;
no strict 'refs';

#use PerlIO::gzip;
my $usage = "
NOTICE: This script works with paired-end fastq data produced by the tang-version of smart-seq.
It will 1)add UMI to sequence_name (line 1 of 4-line fastq) of reads in R1 and 2)split R1 into 
several OUT_PREFIXn.fq.gz files according to the barcode of reads in R2.
AUTHOR: Based on multiple scripts created by Ji Dong, Zongcheng Li shorten the codes into this 
single script, which may also bring increased performance.

USAGE: perl $0 R1 R2 OFS OUT_PREFIX BARCODE_FILE
All parameters are required!

R1: Read_1.fq.gz (Read_1 file of paired-end fastq files)
R2: Read_2.fq.gz (Read_2 file of paired-end fastq files)
OFS: offset bases from the start of a sequence, must less than (read length - 8)
BARCODE_FILE:
BC1
BC2
...

that means 
BC1: BarCode 1, such as GGTCTTAT and so on(coercive: length of barcode must be 8!)
BC2: ditto
...: BCx and so on

OUT_PREFIX: for output files, that will be OUT_PREFIX1.UMI.fq.gz or OUT_PREFIXn.UMI.fq.gz
            including reads extracted from R1 corresponding to BCn 
            and OUT_PREFIX.others.fq.gz (other barcodes)";
die $usage unless @ARGV ==5;
my (@barcodes,$i,$out,$n,%hash);
open IN0,"$ARGV[4]" or die $!;
@barcodes = <IN0>;
$n = @barcodes;

open IN1,"gzip -dc $ARGV[0] |" or die $!;
open IN2,"gzip -dc $ARGV[1] |" or die $!;
for($i=1;$i<=$n;$i++){
         chomp($barcodes[$i - 1]);
	 chomp($barcodes[$i - 1]);
         $hash{$barcodes[$i - 1]} = $i;         
         $out = "OUT".$i;#$out bu gu ding
         `mkdir -p $ARGV[3]$i` unless(-d "$ARGV[3]$i");
         open $out,"| gzip -c > '$ARGV[3]$i'.UMI.fq.gz" or die $!;
         }
`mkdir -p $ARGV[3]"_others"` unless(-d "'$ARGV[3]'_others");
open OUT0,"| gzip -c > '$ARGV[3]'.others.UMI.fq.gz" or die $!;

while (1) {
    #-get the reads and corresponding information in each 4 lines
    
    my $read1_1 = <IN1>;
    my $read1_2 = <IN1>;
    my $read1_3 = <IN1>;
    my $read1_4 = <IN1>;
    
    my $read2_1 = <IN2>;
    my $read2_2 = <IN2>;
    my $read2_3 = <IN2>;
    my $read2_4 = <IN2>;
    
    #check the end of the file
    last unless (defined($read1_1) and defined($read2_1));
    chomp ($read1_1,$read1_2,$read1_3,$read1_4,$read2_1,$read2_2,$read2_3,$read2_4);
    my $seq_UMI = substr($read2_2,8,8);
    
	$read1_1 = "\@$seq_UMI\_$read1_1";
	if (not exists $hash{substr($read2_2,$ARGV[2],8)}){
	#	print OUT0 "$read2_1\n$read2_2\n$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
		print OUT0 "$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
	}else{
	    $out = "OUT".$hash{substr($read2_2,$ARGV[2],8)};
	#    print $out "$read2_1\n$read2_2\n$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
	    print $out "$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
	}
}

close IN1;
close IN2;
close OUT0;
`mv '$ARGV[3]'.others.UMI.fq.gz '$ARGV[3]'_others/`;
for($i=1;$i<=$n;$i++){
         $out = "OUT".$i;
         close $out;
		 `mv '$ARGV[3]$i'.UMI.fq.gz $ARGV[3]$i/`;
        }
