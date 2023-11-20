#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_38,singularity --input assets/test.csv