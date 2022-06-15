#!/usr/bin/perl
use strict;
use Getopt::Long;

Getopt::Long::Configure("bundling");

my $min_len = 16;
my $input = '';
my $output_file ='';
my $adpremove_fastq = 0;
my $bad_read_fastq = 0;
my $short_read_fastq = 0;

GetOptions (
	"input|i=s" => \$input,
	"minimum_read_length|m:i" => \$min_len,
	"output_prefix|o:s" => \$output_file,
	"output_adaptor_removed_fastq|r" => \$adpremove_fastq,
	"retain_bad_read|b" => \$bad_read_fastq,
	"output_short_read_fastq|s" => \$short_read_fastq
);

($output_file = $input) =~ s/.cutadapt.info.txt//;
$output_file =~ s/^.+\///;
unless(open(LOG, ">$output_file.process.log")){die;}

print "Out prefix set to: $output_file\n"; #exit;
print LOG "Out prefix set to: $output_file\n"; #exit;

## Open up the input file, the info-file output from cutadapt
unless(open(IN, "$input")){die;}

## Set the variable for minimal read length to filter out. The default is 12

if ($adpremove_fastq == 1) {
	unless(open(ADPREMOVE, ">$output_file.adaptor_removed.fastq")){die;}
}

if ($short_read_fastq == 1 ) {
	unless(open(SHORT, ">$output_file.adaptor_removed.short_read.fastq")){die;}
}


if ($bad_read_fastq == 1 ) {
	unless(open(BADREAD, ">$output_file.discarded.fastq")){die;}
}

unless(open(DEDUP, ">$output_file.umi_removed.fastq")){die;}



#my @read_array;
my $array_count = 0;

my @a_array; my $a_array_count=0;
my @c_array; my $c_array_count=0;
my @g_array; my $g_array_count=0;
my @t_array; my $t_array_count=0;

my $bad_read_count = 0;
my $short_umi_count = 0;
my $short_read_count = 0;
my $no_adaptor_count = 0;
my $n_umi_count = 0;
my $n_read_count = 0;
my $read_count = 0;
while (<IN>){
	$read_count++;
	if ( $read_count % 250000 == 0) {
		print "Processed $read_count reads\n";
	}
	chomp $_;
	my @line = split(/\t/, $_);	
	$line[0] = "\@$line[0]";
	my $umiplus = $line[6];
	my $umi = '';
	my $read = $line[4];
	my $read_len = length ($read);

	if ( length ($umiplus) > 12){
		$umi = substr ($umiplus, 0, 12);
	}

	if ( length ($umiplus) >= 12 && $line[1] != -1 && $umi !~ m/N/ && $read !~ m/N/ ){
		if ( $read_len >= $min_len) {
			#$read_array[$array_count][0] = $line[0];
			#$read_array[$array_count][1] = $line[4];
			#$read_array[$array_count][2] = $line[8];
			#$read_array[$array_count][3] = $umi;
			$array_count++;

			if ( $umi =~ m/^A/ ){
				$a_array[$a_array_count][0] = $line[0];
				$a_array[$a_array_count][1] = $line[4];
				$a_array[$a_array_count][2] = $line[8];
				$a_array[$a_array_count][3] = $umi;
				$a_array_count++;
			} elsif ($umi =~ m/^C/ ){
				$c_array[$c_array_count][0] = $line[0];
				$c_array[$c_array_count][1] = $line[4];
				$c_array[$c_array_count][2] = $line[8];
				$c_array[$c_array_count][3] = $umi;
				$c_array_count++;
			} elsif ($umi =~ m/^G/ ){
				$g_array[$g_array_count][0] = $line[0];
				$g_array[$g_array_count][1] = $line[4];
				$g_array[$g_array_count][2] = $line[8];
				$g_array[$g_array_count][3] = $umi;
				$g_array_count++;
			} else {
				$t_array[$t_array_count][0] = $line[0];
				$t_array[$t_array_count][1] = $line[4];
				$t_array[$t_array_count][2] = $line[8];
				$t_array[$t_array_count][3] = $umi;
				$t_array_count++;
			}

		
			if ( $adpremove_fastq == 1 ) {
				print ADPREMOVE "$line[0]\n$line[4]\n\+\n$line[8]\n";
			}
		} else {
			print SHORT "$line[0]\n$line[4]\n\+\n$line[8]\n";
			$short_read_count++;
			$bad_read_count++;
		}

	} else {


		if ($line[1] == -1) {
			print BADREAD "$line[0]\n$line[2]\n\+\n$line[3]\n";
			$no_adaptor_count++;
		} elsif ($umi =~ m/N/){
			print BADREAD "$line[0]\n$line[4]\n\+\n$line[8]\n";			
			$n_umi_count++;
		} elsif ( length ($umiplus) < 12 ) {
			print BADREAD "$line[0]\n$line[4]\n\+\n$line[8]\n";
			$short_umi_count++;
		} else {
			print BADREAD "$line[0]\n$line[4]\n\+\n$line[8]\n";
			$n_read_count++;
		}

		$bad_read_count++;
	}	
}

my $sum_array_count = $a_array_count + $c_array_count + $g_array_count + $t_array_count;

print "Out of total $read_count reads\n  Good read: $array_count\n  Discard read : $bad_read_count\n\tno adaptor = $no_adaptor_count\n\tN in UMI = $n_umi_count\n\tN in read = $n_read_count\n\tshort_umi = $short_umi_count\n\tshort_read (< $min_len bases) = $short_read_count\n";
print LOG "Out of total $read_count reads\n  Good read: $array_count\n  Discard read : $bad_read_count\n\tno adaptor = $no_adaptor_count\n\tN in UMI = $n_umi_count\n\tN in read = $n_read_count\n\tshort_umi = $short_umi_count\n\tshort_read (< $min_len bases) = $short_read_count\n";

close IN;
close BADREAD;

print "sorting UMI\n";

print "Sorting sorted UMI started with A\n";
my @sorted_a_array = sort { $a->[3] cmp $b->[3] || $a->[1] cmp $b->[1] } @a_array;
print "Sorting sorted UMI started with C\n";
my @sorted_c_array = sort { $a->[3] cmp $b->[3] || $a->[1] cmp $b->[1] } @c_array;
print "Sorting sorted UMI started with G\n";
my @sorted_g_array = sort { $a->[3] cmp $b->[3] || $a->[1] cmp $b->[1] } @g_array;
print "Sorting sorted UMI started with T\n";
my @sorted_t_array = sort { $a->[3] cmp $b->[3] || $a->[1] cmp $b->[1] } @t_array;

my @sorted_read_array;
my $sorted_read_array_count=0;

print "Adding $a_array_count sorted UMI started with A\n";
for (my $i = 0; $i < $a_array_count; $i++) {
	for (my $j=0 ; $j < 4 ; $j++){
		$sorted_read_array[$sorted_read_array_count][$j] = $sorted_a_array[$i][$j];
	}
	$sorted_read_array_count++;
}

print "Adding $c_array_count sorted UMI started with C\n";

for (my $i = 0; $i < $c_array_count; $i++) {
	for (my $j=0 ; $j < 4 ; $j++){
		$sorted_read_array[$sorted_read_array_count][$j] = $sorted_c_array[$i][$j];
	}
	$sorted_read_array_count++;
}

print "Adding $g_array_count sorted UMI started with G\n";

for (my $i = 0; $i < $g_array_count; $i++) {
	for (my $j=0 ; $j < 4 ; $j++){
		$sorted_read_array[$sorted_read_array_count][$j] = $sorted_g_array[$i][$j];
	}
	$sorted_read_array_count++;
}

print "Adding $t_array_count sorted UMI started with T\n";

for (my $i = 0; $i < $t_array_count; $i++) {
	for (my $j=0 ; $j < 4 ; $j++){
		$sorted_read_array[$sorted_read_array_count][$j] = $sorted_t_array[$i][$j];
	}
	$sorted_read_array_count++;
}

#for (my $i = 0 ; $i < $sorted_read_array_count; $i++){
#	print "$sorted_read_array[$i][3]\n";
	#print "$sorted_read_array[$i][0]\n$sorted_read_array[$i][1]\n$sorted_read_array[$i][2]\n$sorted_read_array[$i][3]\n";
#}
#exit;

print "removing duplicated UMI from $sorted_read_array_count\n";

my @dedup_read_array;
my $dedup_read_array_count=0;
my $duplicated_umi_count=0;
 
## Write out the first record from the sorted read array to dedup read array
for (my $i = 0; $i < 4 ; $i++){
	$dedup_read_array[0][$i] = $sorted_read_array[0][$i];
}
$dedup_read_array_count++;

for (my $i = 1; $i < $sorted_read_array_count; $i++){
	

	if ( $sorted_read_array_count % 250000 == 0) {
		print "Processed $i read, and $dedup_read_array_count sorted reads\n";
	}

	my $curr_umi = $sorted_read_array[$i][3];
	my $prev_umi = $sorted_read_array[$i-1][3];
	my $curr_seq = $sorted_read_array[$i][1];
	my $prev_seq = $sorted_read_array[$i-1][1];	
	my $curr_seq_length = length ($sorted_read_array[$i][1]);
	my $prev_seq_length = length ($sorted_read_array[$i-1][1]);

	#print "$dedup_read_array_count\t$curr_umi\t$curr_seq\t$curr_seq_length\n\n";

	if ( $curr_umi eq $prev_umi && $curr_seq eq $prev_seq ) {

		#print "$prev_umi\t$prev_seq\t$prev_seq_length\n$curr_umi\t$curr_seq\t$curr_seq_length\n\n";
		$duplicated_umi_count++;

	} else {

		#print "$prev_umi\t$prev_seq\t$prev_seq_length\n$curr_umi\t$curr_seq\t$curr_seq_length\n\n";
		$dedup_read_array[$dedup_read_array_count][0] = $sorted_read_array[$i][0];
		$dedup_read_array[$dedup_read_array_count][1] = $sorted_read_array[$i][1];
		$dedup_read_array[$dedup_read_array_count][2] = $sorted_read_array[$i][2];
		$dedup_read_array_count++;



	}
}

print "Number of unique sequence / umi combination: $dedup_read_array_count\n";
print "Number of duplicate found: $duplicated_umi_count\n";


print LOG "Number of unique sequence / umi combination: $dedup_read_array_count\n";
print LOG "Number of duplicate found: $duplicated_umi_count\n";

for (my $i = 0 ; $i < $dedup_read_array_count ; $i++){
	print DEDUP "$dedup_read_array[$i][0]\n";
	print DEDUP "$dedup_read_array[$i][1]\n";
	print DEDUP '+'."\n";
	print DEDUP "$dedup_read_array[$i][2]\n";
}

close DEDUP;
close ADPREMOVE;


exit;


