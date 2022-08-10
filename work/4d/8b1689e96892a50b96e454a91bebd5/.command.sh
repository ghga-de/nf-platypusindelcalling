#!/bin/bash -euo pipefail
check_samplesheet.py \
    samplesheet_cluster.csv \
    samplesheet.valid.csv

cat <<-END_VERSIONS > versions.yml
"NF_PLATYPUSINDELCALLING:PLATYPUSINDELCALLING:INPUT_CHECK:SAMPLESHEET_CHECK":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
