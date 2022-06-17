#!/usr/bin/env python


import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/benchmark samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,tumor,control,reference
    sample1,../data_tests/BAM/NA06984_T.bam,../data_test/BAM/NA06984_N.bam,../data_tests/REF/17.fasta
    sample2,../data_tests/BAM/NA06985_T.bam, , ,

    For an example see:
    https://raw.githubusercontent.com/christopher-hakkaart/testdata/test_benchmark.csv
    """

    sample_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 4
        HEADER = ["sample", "tumor", "control", "genome"]

        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check all sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            sample, tumor, control, genome  = lspl[: len(HEADER)]

            ## Check sample name entries
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)


            ## Check genome entries
            if genome:
                if genome.find(' ') != -1:
                    print_error("Genome entry contains spaces!",'Line', line)
                if len(genome.split('.')) > 1:
                    if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                        print_error("Genome entry does not have extension '.fasta', '.fa', '.fasta.gz' or '.fa.gz'!",'Line', line)


            ## Create sample mapping dictionary = { sample: [ tumor, control, genome] }
            sample_info = [tumor, control, genome]

            if sample not in sample_dict:
                sample_dict[sample] = [sample_info]
            else:
                if sample_info in sample_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "tumor", "control", "genome"]) + "\n")
            for sample in sorted(sample_dict.keys()):
                for idx, val in enumerate(sample_dict[sample]):
                    fout.write(",".join(["{}_{}".format(sample, val[0].lower())] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
