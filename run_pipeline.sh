#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_37,singularity --outdir result4 --input assets/samplesheet_test_37.csv