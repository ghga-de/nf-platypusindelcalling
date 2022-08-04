#!/bin/bash

module load nextflow
nextflow run main.nf -profile dkfz_config,singularity --runTinda false