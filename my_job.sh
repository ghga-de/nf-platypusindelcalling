#!/bin/bash
samtools collate -u -O "$INPUT"| \\
samtools fastq -1 "$OUTPUT1" -2 "$OUTPUT2" -0 /dev/null -s /dev/null -n