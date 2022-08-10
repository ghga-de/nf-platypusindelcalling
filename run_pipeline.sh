#!/bin/bash

module load nextflow/21.10.6
nextflow run main.nf -profile dkfz_cluster,singularity --runTinda false