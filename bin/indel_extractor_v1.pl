#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
# By Ivo Buchhalter
# Hard coded script to extract indels in the snv pipeline
# Extracts: somatic indels -> indels_PID_somtatic_indels_conf_?_to_10.vcf
#           somatic coding indels -> indels_PID_somatic_coding_indels_conf_?_to_10.vcf
# optionally: germline coding indels -> indels_PID_germline_coding_indels_conf_?_to_10.vcf


# Modified 2026-08-12 @kubranarci: I/O error detection (issue #43): a read error makes <FH> return undef
# just like EOF, so each read loop is followed by an explicit ->error check, and closes are checked
# (failed zcat/tabix pipes only show up in the exit status reported by close)

use strict;
use warnings;
use IO::Handle;

# Issue #43: a write error on STDOUT (e.g. disk full) is only reported when STDOUT is closed
END { CORE::close(STDOUT) or do { print STDERR "Error writing to STDOUT: $!\n"; $? ||= 1 } }
use Getopt::Long;

my $infile;
my $minconf=8;
my $region = 0;
my $somout;
my $funcout;
my $ncout;
my $germlineout;
my $bgzip = "bgzip";
my $tabix = "tabix";
GetOptions (	"infile=s"	 		=> \$infile,		# vcf file, can be bgzipped
				"minconf=i"			=> \$minconf,		# minimum confidence score
				"region=s"			=> \$region,		# file with the region when only a region should be extracted
				"somout=s"			=> \$somout,		# Outputfile for somatic indels
				"funcout=s"			=> \$funcout,		# Outputfile for functional somatic indels
				"ncout=s"			=> \$ncout,			# Outputfile for non coding RNAS
				"germlineout=s"		=> \$germlineout,	# germline coding indels
				"bgzip=s"			=> \$bgzip,
				"tabix=s"			=> \$tabix
) or die "Could not get the options!\n";

if($region ne "0" && (!-f $region || $infile !~ /\.gz$/ || !-f $infile.".tbi")){die "region-file: $region is not a valid file or, infile: $infile is not zipped or there is no index for the infile $infile.tbi\n";}
if(!-f $infile){die "The provided infile: $infile is not a valid file\n";}
if(-f $infile && $infile =~ /\.gz$/){open(IN, "zcat $infile |") or die "Could not run zcat on $infile: $!\n";}
else{open(IN, "<$infile") or die "Could not open the infile: $infile\n";}

if($region ne "0"){print "Region file provided: $region\n";}

open(SOM, ">$somout") or die "Could not open the file $somout\n";
open(COD, ">$funcout") or die "Could not open the file $funcout\n";
open(NCO, ">$ncout") or die "Could not open the file $ncout\n";
if (defined $germlineout)
{
	open(GER, ">$germlineout") or die "Could not open the file $germlineout\n";
}

my $head;
while(<IN>)
{
	chomp;
	$head=$_;
	last if($_ =~ /^#CHR/);
}
die "Error reading $infile: $!\n" if (\*IN)->error;  # issue #43
# we stop after the header, so zcat may get SIGPIPE (issue #43)
close(IN) or (($? & 127) == 13) or die "Error reading header of $infile: " . ($! || "exit status " . ($? >> 8)) . "\n";

my @head=split("\t", $head);
my %col;


my $i = 0;
while($i < @head)
{
	if($head[$i] eq "EXONIC_CLASSIFICATION"){$col{"EXONIC_CLASSIFICATION"} = $i; print "EXONIC_CLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "ANNOVAR_FUNCTION"){$col{"ANNOVAR_FUNCTION"} = $i; print "ANNOVAR_FUNCTION in column ", $i+1,"\n";}
	if($head[$i] eq "CLASSIFICATION"){$col{"CLASSIFICATION"} = $i; print "CLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "CONFIDENCE"){$col{"CONFIDENCE"} = $i; print "CONFIDENCE in column ", $i+1,"\n";}
	if($head[$i] eq "REGION_CONFIDENCE"){$col{"REGION_CONFIDENCE"} = $i; print "REGION_CONFIDENCE in column ", $i+1,"\n";}
	$i++;
}

print SOM $head, "\n";
print COD $head, "\n";
print NCO $head, "\n";
if (defined $germlineout)
{
	print GER $head, "\n";
}

if($region ne "0"){open(IN, "$tabix $infile -B $region |") or die "Could not open the file with tabix and regions\n";}
elsif($infile =~ /\.gz$/){open(IN, "zcat $infile |") or die "Could not run zcat on $infile: $!\n";}
else{open(IN, "<$infile") or die "Could not open the infile: $infile\n";}

while(<IN>)
{
	chomp;
	next if($_ =~ /^#/);
	my @line = split("\t", $_);
	if (defined $germlineout && $line[$col{"REGION_CONFIDENCE"}] >= $minconf && ($line[$col{"CLASSIFICATION"}] eq "germline" || $line[$col{"CLASSIFICATION"}] eq "SNP_support_germline") && ($line[$col{"ANNOVAR_FUNCTION"}] eq "exonic" || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/) && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA_exonic/)
	{
		print GER $_, "\n";
	}
	next if($line[$col{"CONFIDENCE"}] < $minconf);
	if($line[$col{"CLASSIFICATION"}] eq "somatic"){print SOM $_, "\n";}
	if($line[$col{"CLASSIFICATION"}] eq "somatic" && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA/ && ($line[$col{"ANNOVAR_FUNCTION"}] =~ /exonic/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/)){print COD $_, "\n";}
	if($line[$col{"CLASSIFICATION"}] eq "somatic" && ($line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_exonic/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_splicing/)){print NCO $_, "\n";}
}
die "Error reading $infile: $!\n" if (\*IN)->error;  # issue #43

close(IN) or die "Error reading $infile: " . ($! || "exit status " . ($? >> 8)) . "\n";  # issue #43
close(SOM) or die "Error writing $somout: $!\n";
close(COD) or die "Error writing $funcout: $!\n";
close(NCO) or die "Error writing $ncout: $!\n";
if (defined $germlineout)
{
	close(GER) or die "Error writing $germlineout: $!\n";
}

