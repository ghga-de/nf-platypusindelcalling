# nf-core/platypusindelcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - 28.05.2024

### `Added`

- assets/config/convertToStdVCF.json and bin/convertToStdVCF.py 
    - Option to output VCF files (all) in standard format (4.2) is added.  Also, TSV formatted confidence annotated/filtrated files are being converted into standard vcf.

- minor changes: 
    - output names
### `Fixed`

- Contig processing is only available for hg38 reference. ALT and/or HLA contigs can be given external in a file. 
    - Automatic generation of HLA/ALT contigs is now possible through tumor BAM instead of fasta.

- conda links in nf-core modules is fixed. 
    - NOTE: Conda enviroments are not available for the pipeline. Holding conda environment.yml links the same in default creates error even when enable_conda is false.  

- bin/confidenceAnnotation_Indels.py
    - flag parsing is generic now. 
    - updated to latest version in https://github.com/DKFZ-ODCF/IndelCallingWorkflow/tree/hg38

- conf/modules.config platypus arguments is fixed (now it is same as dkfz)

### `Dependencies`

### `Deprecated`

- FREQ based filtering is removed from bin/confidenceAnnotation_Indels.py. 

## v1.0dev 

Initial release of nf-core/platypusindelcalling, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
