#!/bin/bash
module load nextflow/23.10.1
nextflow run main.nf -profile dkfz_cluster_38,singularity --input assets/samplesheet_test_38.csv