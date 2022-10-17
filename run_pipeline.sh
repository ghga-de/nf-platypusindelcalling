#!/bin/bash

nextflow run main.nf -profile dkfz_cluster,singularity --outdir results2 -resume --runTinda false