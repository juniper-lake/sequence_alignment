#!/usr/bin/perl

# Name: Juniper Lake
# Class: CISC 636 - Computational Biology and Bioinformatics
# Homework 1, Problem 3: Implement Needleman-Wunsch algorithm for global alignment of DNA sequences

# How to use program: perl lakej_align.txt [-o 1] <input_file>
# where [-o 1] is option to output the dynamic programming table
# Note: it's best to print dynamic programing tables to a file with long sequences


use warnings;
use strict;
use autodie;
use Data::Dumper;

#----------------------------------------------------------#
# store command line args and run main function            #
#----------------------------------------------------------#

my $gap = -8;
my $printoption = 0;  # will print dynamic programming matrix when equal to 1
my $input_file;

# store command line arguments
while (@ARGV)
{
	$_ = shift;
	if (/^-o/) {$_ = shift; if ($_ == 1){$printoption = 1; next;}}
	else {$input_file = $_; next;}
}

# run main function
main();


#-----main-------------------------------------------------#
# compare every seq in input file against all others       #
#----------------------------------------------------------#
sub main
{
	my @seqs;
	my @seqdictionary;
	my $seqinfo = '';
	open(my $input_fh, '<', $input_file);  # read in sequences
	while (<$input_fh>)
	{
		chomp;
		if (/^>/) { $seqinfo = $_; next; } # store sequence descriptions
		elsif (/^\s*$/) {next;} # skip blank lines
		else {push @seqs, $_; push @seqdictionary, [$seqinfo, $_];} # store sequences and descriptions
	}
	# my @seqs = <$input_fh>;
	my $numseqs = scalar(@seqs);
	my $bestx = '';
	my $bestxinfo = '';
	my $bestyinfo = '';
	my $besty = '';
	my @bestalignment = (0,0,0);
	my $bestdpmatrix;
	my $besttraceback;
	for (my $i = 0; $i <= ($numseqs-1); $i++)  # compare every sequence with every other without redundancy
	{
		for (my $j = 0; $j <= ($numseqs-1); $j++)
		{
			if ($j > $i)
			{
				my $x = $seqs[$i];
				my $y = $seqs[$j];
				$x =~ s/\s//gi;  #remove whitespace
				$y =~ s/\s//gi;  #remove whitespace
				my $scoretable_ref = scoring_scheme();
				(my $dpmatrix_ref, my $tracebackmatrix_ref) = build_matrix($x, $y, $scoretable_ref);
				my @alignment = get_alignment($dpmatrix_ref, $tracebackmatrix_ref, $x, $y);
				if ($alignment[2] >= $bestalignment[2])  # store the best alignment and replace if any with higher score come along
				{
					$bestx = $seqs[$i];
					$besty = $seqs[$j];
					@bestalignment = @alignment;
					$bestdpmatrix = $dpmatrix_ref;
					$besttraceback = $tracebackmatrix_ref;
				}
			}

		}
	}
	# print best alignment info after comparing all pairwise sequence alignments
	for (my $i=0; $i <= ($numseqs-1); $i++)
	{
		if ($seqdictionary[$i][1] eq $bestx) {$bestxinfo = $seqdictionary[$i][0];} # find descriptions matching best alignment
		if ($seqdictionary[$i][1] eq $besty) {$bestyinfo = $seqdictionary[$i][0];}
		}
	print "\nBest Alignment Score:\t$bestalignment[2]\nSequence 1:\t$bestalignment[0]\nSequence 2:\t$bestalignment[1]\nSequence 1 ID:\t$bestxinfo\nSequence 2 ID:\t$bestyinfo\n\n";
	if ($printoption == 1)
	{
		print_matrix($bestdpmatrix, $besttraceback, $bestx, $besty);
	}
}

#-----scoring_scheme---------------------------------------#
# build scoring lookup table for substitutions/matches     #
#----------------------------------------------------------#
sub scoring_scheme
{
	my %scoretable;
	my @bases = ('A', 'C', 'T', 'G');
	for (my $i=0; $i <= 3; $i++)
	{	for (my $j=0; $j <=3; $j++)
		{
			if ($i == $j) {$scoretable{$bases[$i]}{$bases[$j]} = 4;} # identity
			elsif (abs (($i**2) - ($j**2)) == 9 ) {$scoretable{$bases[$i]}{$bases[$j]} = -2;} # A/G transition
			elsif (abs (($i**2) - ($j**2)) == 3 ) {$scoretable{$bases[$i]}{$bases[$j]} = -2;} # C/T transition
			else {$scoretable{$bases[$i]}{$bases[$j]} = -3;} # transversion
		}
	}
	# print Dumper(\%scoretable);
	return \%scoretable;
}


#-----build_matrix---------------------------------------#
# build matrix of scores and traceback matrix              #
#----------------------------------------------------------#
sub build_matrix
{
	my $x = $_[0];
	my $y = $_[1];
	my @xarr = split //, $x;
	my @yarr = split //, $y;
	my %scoretable = %{$_[2]};
	my $xlength = scalar(@xarr);
	my $ylength = scalar(@yarr);
	# initialize the dynamic programming matrix with top and left values
	my @dpmatrix;
	$dpmatrix[0][0] = 0;
	for (my $i=1; $i <= ($ylength); $i++) {$dpmatrix[$i][0] = $gap * $i;}
	for (my $j=1; $j <= ($xlength); $j++) {$dpmatrix[0][$j] = $gap * $j;}
	# print Dumper(@dpmatrix);
	# initialize the traceback matrix with top and left values
	my @tracebackmatrix;
	$tracebackmatrix[0][0] = 'X';
	for (my $i=1; $i <= ($ylength); $i++) {$tracebackmatrix[$i][0] = '|';}  # all vertical moves on left
	for (my $j=1; $j <= ($xlength); $j++) {$tracebackmatrix[0][$j] = '-';} # all horizontal moves on top
	# fill the other values in the matrix
	for (my $i=1; $i <= ($ylength); $i++)
		{for (my $j=1; $j <= ($xlength); $j++)
			{
				my $h = $dpmatrix[$i][$j-1] + $gap;  # horizontal move
				my $v = $dpmatrix[$i-1][$j] + $gap;  # vertical move
				my $d = $dpmatrix[$i-1][$j-1] + $scoretable{$yarr[$i-1]}{$xarr[$j-1]};  # diagonal move
				my @tmp = ([$v,"|"], [$h,"-"], [$d,"\\"]); 
				# print Dumper(\@tmp);
				my @tmpsorted = sort { $b->[0] <=> $a->[0] } @tmp;  # sort arrays based on scores
				# print Dumper(\@tmpsorted);
				$dpmatrix[$i][$j] = $tmpsorted[0][0];
				$tracebackmatrix[$i][$j] = $tmpsorted[0][1]
			}
		}
		return (\@dpmatrix, \@tracebackmatrix);
}


#-----print_matrix-----------------------------------------#
# print matrices of scores and traceback                   #
#----------------------------------------------------------#
# print the dynamic programming matrix
sub print_matrix
{
	my @dpmatrix = @{$_[0]};
	my @tracebackmatrix = @{$_[1]};
	my $x = $_[2];
	my $y = $_[3];
	my @xarr = split //, $x;
	my @yarr = split //, $y;
	my $i = 0;
	print "Dynamic Programming Score Matrix\n";
	print "\t-\t";
	foreach my $base (@xarr) {print $base,"\t";} # print sequence 1 on top line
	# print "\n";
	foreach my $row (@dpmatrix) {
		if ($i > 0) {print $yarr[$i-1],"\t";}
		else {print "-\t";}
		$i++;
		foreach my $element (@$row) {
			print $element, "\t";

		}
		print "\n";
	}
	print "\n";
	# print the dynamic programming matrix
	$i = 0;
	print "Dynamic Programming Traceback Matrix\n";
	print "\t-\t";
	foreach my $base (@xarr) {print $base,"\t";}
	# print "\n";
	foreach my $row (@tracebackmatrix) {
		if ($i > 0) {print $yarr[$i-1],"\t";}
		else {print "-\t";}
		$i++;
		foreach my $element (@$row) {
			print $element, "\t";

		}
		print "\n";
	}
}


#-----get_alignment----------------------------------------#
# use dynamic programming matrices to visualize alignment  #
#----------------------------------------------------------#
# @dpmatrix = @{$dpmatrix_ref};
sub get_alignment
{
	my @dpmatrix = @{$_[0]};
	my @tracebackmatrix = @{$_[1]};
	my $x = $_[2];
	my $y = $_[3];
	my @xarr = split //, $x;
	my @yarr = split //, $y;
	# print Dumper(\@tracebackmatrix);
	my $score = $dpmatrix[-1][-1];
	my $xalign = '';
	my $yalign = '';
	my $i = -1;
	my $j = -1;
	# print Dumper($tracebackmatrix[$i][$j]);
	while ($tracebackmatrix[$i][$j] ne 'X') 
	{
		if ($tracebackmatrix[$i][$j] eq '\\') {$xalign=$xarr[$j].$xalign; $yalign=$yarr[$i].$yalign; $i--; $j--;}
		elsif ($tracebackmatrix[$i][$j] eq '|') {$xalign='-'.$xalign; $yalign=$yarr[$i].$yalign; $i--;}
		elsif ($tracebackmatrix[$i][$j] eq '-') {$xalign=$xarr[$j].$xalign; $yalign='-'.$yalign; $j--;}
	}
	return ($xalign, $yalign, $score);
	# print $xalign, "\n";
	# print $yalign, "\n";
}

