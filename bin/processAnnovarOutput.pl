#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#
# Modified: 2026-08-12 @kubranarci
# Added: Explicit I/O error handling for nested while loops reading from two input files
#        to ensure consistent and verifiable data processing
#

use strict;
use warnings;
use v5.10;

(@ARGV >= 2) || die "Usage: processAnnovarOutput.pl <variant_function file> <exonic_variant_function file>";
my ($file1, $file2) = @ARGV;
my (@f2_line, @f1_line);
my $i;
my $nr = 0;
#  my $f2lnr; #file 2 line number
my @f2anno;
open(F1, $file1) || die "Could not open variant_function file: $!\n";
open(F2, $file2) || die "Could not open exonic_variant_function file: $!\n";
while (1) {
    my $line = <F2>;
    if (!defined $line) {
        last;  # reached end of F2
    }
    chomp($line);
    @f2_line = split(/\t/, $line);
    #    $f2lnr = substr($line[0],4);
    @f2anno = ($f2_line[1], $f2_line[2]);
    while (1) {
        my $f1_line_str = <F1>;
        if (!defined $f1_line_str) {
            die "Error reading from variant_function file: $!\n" if $!;
            # End of F1 reached, but we still have lines in F2
            say join "\t", $f2_line[3], '.', '.', '.', '.', '.', '.';
            last;
        }
        chomp($f1_line_str);
        #      $nr++;
        @f1_line = split(/\t/, $f1_line_str);
        if ($f2_line[3] eq $f1_line[2] && $f2_line[4] == $f1_line[3] && $f2_line[5] == $f1_line[4]) {
            say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], @f2anno;
            last;
        } else {
            say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
        }
    }
}
while (1) {
    my $line = <F1>;
    if (!defined $line) {
        last;  # reached end of F1
    }
    # print out remaining lines from file 1 after file 2 has ended
    chomp($line);
    @f1_line = split(/\t/, $line);
    say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
}
close(F1);
close(F2);
