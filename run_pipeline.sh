#!/bin/bash

module load anaconda3/2021.05
conda activate /dkfz/cluster/virtualenvs/paramasi/nf
nextflow run main.nf -profile dkfz_cluster,singularity -resume
