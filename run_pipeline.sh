#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile test,singularity --input testdata/samplesheet_test.csv