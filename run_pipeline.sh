#!/bin/bash

nextflow run main.nf -profile dkfz_cluster,singularity --outdir result_test -resume