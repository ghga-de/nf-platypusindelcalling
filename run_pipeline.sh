#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_38,singularity --input testdata/samplesheet_test.csv -resume