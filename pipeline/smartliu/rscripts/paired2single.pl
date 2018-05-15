#!/usr/bin/env perl
use strict;
use warnings;
use JSON;
use utf8;
use Encode;
no strict 'refs';
# perl paired2single.pl ivs_test_1.fq.gz ivs_test_2.fq.gz 96-8bp-barcode 1 ivs_test.single.fq.gz ivs_test.others.fq.gz TGGTATCAACGCAGAGTACAT AAAAAAAAAAAAAAA 40

my $usage = "perl $0 <fq1s> <fq2s> <barcodes> <mismatch> <OUT1> <OUT2> <TSO> <PolyA> <min_length> <json>
<IN1>:fastq.gz
<OUT>:fastq.gz
";
die $usage unless $#ARGV == 9;
my @in1 = split(/,/, $ARGV[0]);
my @in2 = split(/,/, $ARGV[1]);

if (@in1 ne @in2){
	print 'R1/R2 fastq file is not paired!';
	exit();
}

my @barcodes = split(/,/, $ARGV[2]);
my $n = @barcodes;
my $mismatch = $ARGV[3];
my %mis_dict = mismatch_dict(\@barcodes, $mismatch);
#print join(",",values(%mis_dict));
open OUT1,"| gzip -c > $ARGV[4]" or die $!;
open OUT2,"| gzip -c > $ARGV[5]" or die $!;
my $tso = $ARGV[6];
my $tso_n = length($tso);
my $polya = $ARGV[7];
my $min_len = $ARGV[8];
open OUT0,"> $ARGV[9]" or die $!;

my %bbcount = init_bbcount(\@barcodes, $mismatch);


for my $fq (0..$#in1){
	open IN1,"gzip -dc $in1[$fq] |" or die $!;
	open IN2,"gzip -dc $in2[$fq] |" or die $!;
	my $isbar_flag = -1;
	my $ind_flag = -1;
	my $ind_tso = -1;
	my $ind_polya = -1;
	my @bb;
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
		
		$isbar_flag = 0;
		chomp ($read1_1,$read1_2,$read1_3,$read1_4,$read2_1,$read2_2,$read2_3,$read2_4);
		my $tag = substr($read2_2,0,8);
		if (exists $mis_dict{$tag}){
			@bb = @{$mis_dict{$tag}};
			$bbcount{$bb[1]}[0][$bb[0]] += 1;
			$isbar_flag = 1;
		}else{
			$bbcount{'unmatched'}[0] +=1;
		}
		
		$ind_tso = 0;
		$ind_polya = length($read1_2);
		$ind_flag = 0;
		if ($read1_2 =~ m/$tso/g){
			$ind_tso = rindex($read1_2, $tso) + $tso_n + 3;
			$ind_flag = 1;
		}
		if($read1_2 =~ m/$polya/g){
			$ind_polya = index($read1_2, $polya);
			$ind_flag = 1;
		}
		
		if ($ind_flag eq 1){
			if(($ind_polya - $ind_tso) >= $min_len){
			  $read1_2 = substr($read1_2, $ind_tso, $ind_polya - $ind_tso);
				$read1_4 = substr($read1_4, $ind_tso, $ind_polya - $ind_tso);
			}else{
			  next;
			}
		}
		
		if ($isbar_flag eq 1){
			$read1_1 = '@'."$bb[1]".substr($read2_2, 8, 8)."_"."$read1_1";
			$bbcount{$bb[1]}[1][$bb[0]] += 1;
			print OUT1 "$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
		}else{
			$read1_1 = '@'.substr($read2_2, 0, 16)."_"."$read1_1";
			$bbcount{'unmatched'}[1] += 1;
			print OUT2 "$read1_1\n$read1_2\n$read1_3\n$read1_4\n";
		}
	}
	close IN1;
	close IN2;
}
close OUT1;
close OUT2;
my $json = encode_json \%bbcount;
#print "$json";
print OUT0 "$json";
close OUT0;

############Functions############
sub hamming2 {
	my ($s1, $s2) = @_;
	my $error = 0;
	if (length($s1) ne length($s2)){
		return 100;
	}else{
		my $i = 0;
		while ($i <= length($s1) - 1){
			if (substr($s1,$i,1) ne substr($s2,$i,1)){
				$error += 1;
			}
			$i += 1;
		}
		return $error;
	}
}

sub mismatch_dict {
	my ($barcode, $mis)=@_;
	my @base=('A','T','G','C','N');
	my %barcodes_mis_dict;
	my $bar_mis;
	foreach my $bar (@$barcode){
		if ($mis eq 0){
			$barcodes_mis_dict{$bar} = $bar;
		}elsif ($mis eq 1){
			for (my $i=0; $i<length($bar); $i++){
				foreach my $b (@base){
				  if($i eq 0){
  					$bar_mis = "$b".substr($bar, $i+1,length($bar)-$i-1);
				  }elsif($i eq 7){
  					$bar_mis = substr($bar,0,$i)."$b";
				  }else{
  					$bar_mis = substr($bar,0,$i)."$b".substr($bar, $i+1,length($bar)-$i-1);
				  }
					$barcodes_mis_dict{$bar_mis} = [hamming2($bar_mis, $bar), $bar];
					#print $barcodes_mis_dict{$bar_mis}[1];
				}
			}
		}elsif ($mis eq 2){
			for (my $i=0; $i<length($bar)-1; $i++){
				for (my $j=$i+1; $j<length($bar); $j++){
					foreach my $b1 (@base){
						foreach my $b2 (@base){
						  if($i eq 0){
  				    	$bar_mis = "$b1".substr($bar, $i+1,$j-$i-1)."$b2".substr($bar, $j+1, length($bar)-$j-1);
    				  }elsif($j eq 7){
      					$bar_mis = substr($bar,0,$i)."$b1".substr($bar, $i+1,$j-$i-1)."$b2";
    				  }else{
      					$bar_mis = substr($bar,0,$i)."$b1".substr($bar, $i+1,$j-$i-1)."$b2".substr($bar, $j+1, length($bar)-$j-1);
    				  }
							$barcodes_mis_dict{$bar_mis} = [hamming2($bar_mis, $bar), $bar];
						}
					}
				}
			}
		}else{
			print 'mismatch should be less than 3 when matching 8bp barcodes!';
			exit();
		}
	}
	return %barcodes_mis_dict
}

sub init_bbcount {
	my ($barcode, $mis)=@_;
	my %hash;
	$hash{'unmatched'} = [(0,0)];
	foreach my $bar (@$barcode){
		my (@count1, @count2);
		for my $m (0..$mis){
			push(@count1, 0);
			push(@count2, 0);
		}
		$hash{$bar} = [[@count1], [@count2]];
		#print $hash{$bar}."\n";
		#print $hash{$bar}[1][0];
		#print join(",", $hash{$bar}[1])."\n";
	}
	return %hash
}