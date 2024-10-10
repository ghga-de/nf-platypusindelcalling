#!/bin/bash
module load nextflow/23.10.1
nextflow run main.nf -profile dkfz_cluster_38,singularity --input testdata/samplesheet_test.csv -resume