#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile test2,singularity --outdir seq2c_small --input assets/samplesheet_test.csv