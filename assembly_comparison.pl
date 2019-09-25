#/usr/bin/perl -w
use strict;
use Data::Dumper;

#Parses mashmap.out and gets number of differences between genome assemblies
#Joana Damas
#24 September 2019

# check the number of argument
if ($#ARGV+1 != 3) {
	print STDERR "Couldn't find all arguments!\nUsage: perl parse_mashmap_out_v2.pl <mashmap out> <reference assembly sizes> <mashmap resolution>\n";
	exit(1);
}

my $input = $ARGV[0]; #mashmap.out file
my $sizes = $ARGV[1]; #reference scaffold sizes. Format: scaffoldID\tsize
my $resolution = $ARGV[2]; #resolution used for mashmap, in bp.

my $output = "${input}.diff"; #output file without extension
my $new_res = $resolution - 1;
$resolution = $new_res; 

###############################
## START READING INPUT FILES ##
###############################

my %data;
open (IN, $input) or die "Couldn't open $input!\n";
while (<IN>) {
	chomp;
	my @tmp = split(/\s+/, $_);
	my $tarid = $tmp[0];
	if (exists $data{$tarid}) {
		push(@{$data{$tarid}}, join("*", @tmp));
	}
	else {
		@{$data{$tarid}} = join("*", @tmp);
	}
}
close IN;

my %sizes;
open (S, $sizes) or die "Couldn't open $sizes";
while(<S>){
	chomp;
	my @tmp = split(/\t/, $_);
	my $scaf = $tmp[0];
	my $len = $tmp[1];
	$len =~ s/,//g;
	$sizes{$scaf} = $len;
}
close S;

######################################
##     END READING INPUT FILES      ##
######################################
## START MERGING OVERLAPPING BLOCKS ##
######################################

my %eh_data;
my $cnt = 0;
foreach my $key (keys %data){
	foreach my $val (@{$data{$key}}){
		my @tmp = split(/\*/, $val);
		my $lentar = $tmp[3] - $tmp[2]; 
		my $lenref = $tmp[8] - $tmp[7];
		if ($lenref >= $resolution && $lentar >= $resolution){
			my $or = "+1";
			if ($tmp[4] eq "-"){ $or = "-1"; }
			my $nstart = $tmp[2];
			my $nend = $tmp[3];
			my $nrstart = $tmp[7];
			my $nrend = $tmp[8];

			my $chr1 = $tmp[5];
			my $chr2 = $tmp[0];

			@{$eh_data{$cnt}}=("refID", $chr1, $nrstart, $nrend, $nstart, $nend, $or, "tarID", $chr2, $chr2);
            $cnt++;
		}
	}
}

my @merged = mergeBlocks(%eh_data, $resolution);

####################################
## END MERGING OVERLAPPING BLOCKS ##
####################################
##   START FINGDING DIFFERENCES   ##
####################################

my @letter_matrix = &add_letter(@merged);

#Get number of fragments for each reference scaffold
my %frags;
for my $row (0..$#letter_matrix) {
	my @tmp = @{$letter_matrix[$row]};
	if (defined $tmp[11] && $tmp[11] ne ""){
		if(exists $frags{$tmp[7]}{$tmp[1]}){
			push (@{$frags{$tmp[7]}{$tmp[1]}}, $tmp[11]);
		}
		else{
			@{$frags{$tmp[7]}{$tmp[1]}} = $tmp[11];
		}
	}
}
my @sorted_matrix = sort { $a->[7] cmp $b->[7] || $a->[8] cmp $b->[8] || $a->[4] <=> $b->[4] } @letter_matrix; #sort by compGen, comScaf, compScafStart
my @sorted_mat = sort { $a->[7] cmp $b->[7] || $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @letter_matrix; #sort by compGen, refScaf, refScafStart

my %summary;
my %support;
for my $i (0 .. scalar(@sorted_matrix)-2){
	my @tmp1 = @{$sorted_matrix[$i]};
	my @tmp2 = @{$sorted_matrix[$i+1]};
	my @repeated;
	if ($tmp1[8] eq $tmp2[8]){
		my ($suffix1, $suffix2);
		for my $row (0..$#sorted_matrix) { #get frag id 
			if(@{$sorted_matrix[$row]}[10] eq $tmp1[10]){ $suffix1 = @{$sorted_matrix[$row]}[11]; }
			if(@{$sorted_matrix[$row]}[10] eq $tmp2[10]){ $suffix2 = @{$sorted_matrix[$row]}[11]; }
		}
		my ($class1, $class2);
		if ($tmp1[6] eq "+1"){ $class1 = &check_end($tmp1[1], ${suffix1}, $sizes{$tmp1[1]}, $tmp1[2], $tmp1[3], $tmp1[6], $tmp1[7]); }
		if ($tmp1[6] eq "-1"){ $class1 = &check_start($tmp1[1], ${suffix1}, $tmp1[2], $tmp1[3], $tmp1[6]); }
		if ($tmp2[6] eq "+1"){ $class2 = &check_start($tmp2[1], ${suffix2}, $tmp2[2], $tmp2[3], $tmp2[6]); }
		if ($tmp2[6] eq "-1"){ $class2 = &check_end($tmp2[1], ${suffix2}, $sizes{$tmp2[1]}, $tmp2[2], $tmp2[3], $tmp2[6], $tmp2[7]); }

		if (defined $suffix1){
			push (@repeated, ($tmp1[7], $tmp1[8], $tmp1[4], $tmp1[5], $tmp2[4], $tmp2[5], $tmp1[0], $tmp1[1], ${suffix1}, $sizes{$tmp1[1]}, $tmp1[2], $tmp1[3], $tmp1[6]));
		} else { 
			push (@repeated, ($tmp1[7], $tmp1[8], $tmp1[4], $tmp1[5], $tmp2[4], $tmp2[5], $tmp1[0], $tmp1[1], "", $sizes{$tmp1[1]}, $tmp1[2], $tmp1[3], $tmp1[6]));
		}
		if (defined $suffix2){
			push (@repeated, ($tmp2[1], ${suffix2}, $sizes{$tmp2[1]}, $tmp2[2], $tmp2[3], $tmp2[6], $class1, $class2));
		} else {
			push (@repeated, ($tmp2[1], "", $sizes{$tmp2[1]}, $tmp2[2], $tmp2[3], $tmp2[6], $class1, $class2));
		}
	} else{ next; }

	#print Dumper @repeated;	
	&check_for_repetitions(@repeated);

}

##################################
##   END FINGDING DIFFERENCES   ##
##################################
##      PRINT OUTPUT FILES      ##
##################################

open (OUT, ">$output") or die "Couldn't create $output!\n";
print OUT "RefGen\tRefScaf1ID\tRefScaf1Frag\tRefScaf1Len\tRefScaf1Start\tRefScaf1End\tRefScaf1Or\tRefScaf1Adj\t";
print OUT "RefScaf2ID\tRefScaf2Frag\tRefScaf2Len\tRefScaf2Start\tRefScaf2End\tRefScaf2Or\tRefScaf2Adj\t";
print OUT "CompGen\tCompScafID\tCompScaf1Start\tCompScaf1End\tCompScaf2Start\tCompScaf2End\tDiffType\n";
my @classifications;
for my $key1(keys %summary){
	for my $key2 (keys %{$summary{$key1}}){
		my ($type, $ebrtype);
		my @array = split("\t", @{$summary{$key1}{$key2}}[0]);
		if ($array[7] =~ "middle" || $array[14] =~ "middle") {
			$type = "break";
		} else { $type = "end2end"; }
		if ($array[1] eq $array[8]) {
			$ebrtype = "intra";
		} else { $ebrtype = "inter"; }
		push (@{$summary{$key1}{$key2}}, "${type}_${ebrtype}\n");
		push (@classifications, $type);
		print OUT @{$summary{$key1}{$key2}};
	}
}
close OUT;

my $end2end_cnt = grep { $_ eq 'end2end' } @classifications;
my $incons_cnt = grep { $_ eq 'break' } @classifications;

#GET BROKEN CONTIGS/SCAFFOLDS
my @broken;
open (OUT4, ">${output}.broken") or die "Couldn't create ${output}.broken!\n";
print OUT4 "tarID\ttarScafID\ttarScafBreakStart\ttarScafBreakEnd\trefScaf1ID\trefScaf2ID\n";
for my $i (0 .. $#sorted_mat - 1){
	my @tmp1 = @{$sorted_mat[$i]};
	my @tmp2 = @{$sorted_mat[$i+1]};
	if ($tmp1[1] eq $tmp2[1] && $tmp1[8] ne $tmp2[8]){ #REFCONTIG is the same, but COMPSCAF is different
		print OUT4 $tmp1[7]."\t".$tmp1[1]."\t".$tmp1[3]."\t".$tmp2[2]."\t".$tmp1[8]."\t".$tmp2[8]."\n";
		push (@broken, $tmp1[1]);
	}
}	
close OUT4;

my $no_break = scalar(@broken);
my @uniq_scafs = uniq(@broken);
my $no_scafs = scalar(@uniq_scafs);

open (O, ">${output}.summary") or die "Couldn't create ${output}.summary!\n";
print O "No. end to end joins = $end2end_cnt\n";
print O "No. inconsistent joins = $incons_cnt\n";
print O "No. breaks as free ends\n";
print O "    No. breaks = $no_break\n";
print O "    No. scaffolds = $no_scafs\n";
close O;

#########
## END ##
#########

##SUBROUTINES
sub mergeBlocks {
    my (%data_ref, $resolution)=@_;
    #my %data = %$data_ref;
    my @merged_data;
    my @sorted_keys=();
    for my $key ( sort { $data_ref{$a}[1] cmp $data_ref{$b}[1] || $data_ref{$a}[2] <=> $data_ref{$b}[2] } keys %data_ref ) {
        push(@sorted_keys, $key);  
    }

    my $add;
    my $row = 0;
    START: for ( my $k = 0; $k <= scalar (@sorted_keys); $k += $add){ #for each line in file ordered by ref chr and ref start
        my ($label11, $chr1, $chrS1, $chrE1, $carS1, $carE1, $carOr1, $label21, $carI1) = @{$data_ref{$sorted_keys[$k]}};
        my ($bstart, $bend, $cstart, $cend)=($chrS1, $chrE1, $carS1, $carE1);
        GO: for (my $i = 1; $i <= scalar (@sorted_keys) ; $i ++){
	        if (($k+$i) == scalar(@sorted_keys)) { #if last line of file
                @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
				$cnt++;
                last START;
            }
            my ($label12, $chr2, $chrS2, $chrE2, $carS2, $carE2, $carOr2, $label22, $carI2) = @{$data_ref{$sorted_keys[$k+$i]}};
      
            if ($carOr1 eq "+1" && $chr1 eq $chr2 && $carI1 eq $carI2 && $carOr1 eq $carOr2){ ## IF BOTH BLOCKS + ORIENTATION
                if ( $cend == $carS2 ){ ##IF SUCCESSIVE MERGE
                    $bend = $chrE2;
                    $cend = $carE2;
                    next GO;
                }
                my $gap = $cend - $carS2; # calculate gap in target species (end of 1st block - start of 2nd block)
                if ($gap > 0){ ## NOT SUCCESSIVE -> BREAK (start of 2nd block will be smaller than end of 1st block)
                    if ($carS2 < $cend && $carS2 > $cstart && $carE2 > $cstart){ #except if overlapping at begginig of second block, then merge
                        $bend = $chrE2;
                        $cend = $carE2;
                        next GO;
                    } else{      
                        @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
						$cnt++;
                        $add = $i;
                        next START;
                    }
                }
                if ( $gap < 0 && abs($gap) > $resolution ){ ## MAY BE LACK OF ALIGNMENT -> CHECK FOR BLOCKS IN GAP REGION
                    my $flag = "False";
                    for (my $j = 1; $j < scalar (@sorted_keys) ; $j ++){
                        my @array = ($k .. $k+$i);
                        if ( ! grep (/$j/, @array) ){
                            my ($label13, $chr3, $chrS3, $chrE3, $carS3, $carE3, $carOr3, $label23, $carI3) = @{$data_ref{$sorted_keys[$j]}};
                            if ( $carI2 eq $carI3 && $cend <= $carE3 && $carS3 <= $carS2 ){  $flag = "True"; } # OVERLAP FOUND
                        }             
                    }
                    if ($flag eq "True"){ ## IF OVERLAP FOUND -> BREAK
                        @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
						$cnt++;
                        $add = $i;
                        next START;
                    }else{ ## IF NO OVERLAP -> MERGE
                        $bend = $chrE2;
                        $cend = $carE2;
                        next GO;
                    }
                }
                if ( $gap < 0 && abs($gap) <= $resolution ){ ## IF GAP < RESOLUTION -> MERGE
                    $bend = $chrE2;
                    $cend = $carE2;
                    next GO;
                }
            }
            if ($carOr1 eq "-1" && $chr1 eq $chr2 && $carI1 eq $carI2 && $carOr1 eq $carOr2){ ## IF BOTH BLOCKS INVERTED 
                if ($cstart == $carE2){ ## IF BLOCKS ARE CONSECUTIVE MERGE
                    $bend = $chrE2;
                    $cstart = $carS2;
                    next GO;
                }
                my $gap = $carE2 - $cstart; ## GAP LENGTH
                if ($gap > 0){ ## IF GAP > 0: 2 INDEPENDENT INVERTIONS -> BREAK
                    if($carE2 > $cstart && $cend > $carS2 && $cend > $carE2){ # except if overlapping
                        $bend = $chrE2;
                        $cstart = $carS2;
                        next GO;
                    } else{
                        @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
						$cnt++;
                        $add = $i;
                        next START;
                    }
                }
                if ($gap < 0 && abs($gap) > $resolution){ ## IF GAP < 0: MAY BE LACK OF ALIGNMENT -> CHECK FOR OTHER BLOCKS IN GAP
                    my ($gapStart, $gapEnd)=($carE2, $cstart);
                    my $flag = "False";
                    for (my $j = 1; $j < scalar (@sorted_keys) ; $j ++){
                        my @array = ($k .. $k+$i);
                        if ( ! grep (/$j/, @array) ){
                            my ($label13, $chr3, $chrS3, $chrE3, $carS3, $carE3, $carOr3, $label23, $carI3) = @{$data_ref{$sorted_keys[$j]}};
                            if ( $carI3 eq $carI2 && $gapStart <= $carE3 && $carS3 <= $gapEnd ){ $flag = "True"; } ##OVERLAP FOUND
                        }             
                    }
                    if ($flag eq "True"){ ## IF OVERLAP FOUND -> BREAK
                        @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
						$cnt++;
                        $add = $i;
                        next START;
                    }else{ # OTHERWISE MERGE
                        $bend = $chrE2;
                        $cstart = $carS2;
                        next GO;
                    }
                }
                if ($gap < 0 && abs($gap) <= $resolution ){ ## IF GAP LOWER THAN RESOLUTION -> MERGE
                    $bend = $chrE2;
                    $cstart = $carS2;
                    next GO;
                }  
            }
            if ($chr1 ne $chr2 || $carI1 ne $carI2 || $carOr1 ne $carOr2){
                @{$merged_data[$cnt]} = ($label11, $chr1, $bstart, $bend, $cstart, $cend, $carOr1, $label21, $carI1, $carI1, $cnt);
				$cnt++;
                $add = $i;
                next START;
            }
        }
    }
    #print Dumper @merged_data;
    @merged_data = grep { $_ ne '' } @merged_data;
    return @merged_data;
}

sub add_letter{
	my @contents = @_;
	my @contents_sorted = sort { $a->[7] cmp $b->[7] || $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @contents; #order by compGen, refScaf, refScafStart
	my $prevTag = "";
	my $currentTag = "";
	my $prevSP = "";
	my $currentSP = "";
	my $sum = 0;
	my $hadFirst = '0';
	my $m = $#contents_sorted;
	#print "$m\n";
	my $n = 11; #number of columns

	for my $i (0..$m) {
		$currentTag = $contents_sorted[$i][1];
		$currentSP = $contents_sorted[$i][7];
		#print "$currentTag\n";
		if ($currentTag eq $prevTag and $currentSP eq $prevSP) {
			if ($hadFirst == '1') {
				push (@{$contents_sorted[$i]},chr(ord('a') + $sum));
				$sum++;
			}
			else {
				push (@{$contents_sorted[$i-1]}, chr(ord('a') + $sum));
				$sum++;
				push (@{$contents_sorted[$i]}, chr(ord('a') + $sum));
				$sum++;
				$hadFirst = '1';
			}
		}
		else {
			$hadFirst = '0';
			$sum = 0;
		}
		$prevTag = $currentTag;
		$prevSP = $currentSP;
	}
	return @contents_sorted;
}

sub check_start{
	my $type;
	my ($scaf, $frag, $start, $end, $dir) = @_;
	if ($start == 1 || $start == 0){ $type = "start"; }
	else{
		if (! defined $frag || $frag eq "") { $type = "start"; }
		else{
			if ($frag eq "a"){ $type = "start"; }
			else{ 
				$type = "${start}_middle_start";
				# if ($dir eq "-1" || $dir eq "-"){
				# 	$type = "${start}_middle_start";
				# }
			}
		}
	}
	return $type; 
}

sub check_end{
	my $type;
	my ($scaf, $frag, $len, $start, $end, $dir, $sp) = @_;
	if ( $end == $len || $end == $len-1 ){ $type = "end"; }
	else{
		if (exists $frags{$sp}{$scaf}){
			#print Dumper @{$frags{$scaf}};
			my $noFrag = scalar(@{$frags{$sp}{$scaf}});
			my $lastFrag = chr(ord('a') + $noFrag - 1);
			#print $scaf."\t".$noFrag."\t".chr(ord('a') + $noFrag - 1)."\n";
			if ($frag eq $lastFrag){ 
				$type = "end";
				#print "FRAG ==> $frag\n";
			}
			else{ 
				$type = "${end}_middle_end"; 
				# if ($dir eq "-1" || $dir eq "-"){
				# 	$type = "${end}_middle_end";
				# }	
			}
		}
		else { $type = "end"; }
	}
	return $type;
}

sub check_for_repetitions{
	my @tmp=@_;

	my ($class1, $class2) = ($tmp[19], $tmp[20]);
	#print "$class1\t$class2\n";
	if (exists $summary{$tmp[7].":".$class1}{$tmp[13].":".$class2}) {
		my $line = &create_line(@tmp); 
		push( @{$summary{$tmp[7].":".$class1}{$tmp[13].":".$class2}}, $line); 
		push( @{$support{$tmp[7].":".$class1}{$tmp[13].":".$class2}}, $tmp[0]);
	}
	elsif (exists $summary{$tmp[13].":".$class2}{$tmp[7].":".$class1}) {
		my $line = &create_line_inverse(@tmp); 	
		push( @{$summary{$tmp[13].":".$class2}{$tmp[7].":".$class1}}, $line);
		push( @{$support{$tmp[13].":".$class2}{$tmp[7].":".$class1}}, $tmp[0]);
	}
	else{
		my $line = &create_line(@tmp);
		@{$summary{$tmp[7].":".$class1}{$tmp[13].":".$class2}} = $line;
		@{$support{$tmp[7].":".$class1}{$tmp[13].":".$class2}} = $tmp[0];
	}
} 

sub create_line{
	my (@tmp) = @_;
	my $line;
	my ($suffix1, $suffix2) = ($tmp[8], $tmp[14]);
	if (defined $suffix1){
		$line = "$tmp[6]\t$tmp[7]\t$suffix1\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\t$tmp[19]"; 
	} else{ 
		$line = "$tmp[6]\t$tmp[7]\t\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\t$tmp[19]";
	}
	if (defined $suffix2){
		$line .= "\t$tmp[13]\t$suffix2\t$tmp[15]\t$tmp[16]\t$tmp[17]\t$tmp[18]\t$tmp[20]";
	} else {
		$line .= "\t$tmp[13]\t\t$tmp[15]\t$tmp[16]\t$tmp[17]\t$tmp[18]\t$tmp[20]";
	}
	$line .= "\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t";
	return $line;
}

sub create_line_inverse{
	my (@tmp) = @_;
	my $line;
	my ($suffix1, $suffix2) = ($tmp[14], $tmp[8]);
	if (defined $suffix1){
		$line = "$tmp[6]\t$tmp[13]\t$suffix2\t$tmp[15]\t$tmp[16]\t$tmp[17]\t$tmp[18]\t$tmp[20]";
	} else{ 
		$line = "$tmp[6]\t$tmp[13]\t\t$tmp[15]\t$tmp[16]\t$tmp[17]\t$tmp[18]\t$tmp[20]";
	}
	if (defined $suffix2){
		$line .= "\t$tmp[7]\t$suffix1\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\t$tmp[19]"; 
	} else {
		$line .= "\t$tmp[7]\t\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\t$tmp[19]";
	}
	$line .= "\t$tmp[0]\t$tmp[1]\t$tmp[4]\t$tmp[5]\t$tmp[2]\t$tmp[3]\t";
	return $line;
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}