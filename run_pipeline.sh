#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_38,singularity --outdir result --input assets/samplesheet_test_38.csv