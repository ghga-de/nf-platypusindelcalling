# nf-core/platypusindelcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/platypusindelcalling, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## v2.0.0 - [09.10.2024]

### `Added`

- assets/config/convertToStdVCF.json and bin/convertToStdVCF.py

    - Option to output VCF files (all) in standard format (4.2) is added. Also, TSV formatted confidence annotated/filtrated files are being converted into standard VCF.

- MAFCommon tag is added to the INFO column to mark the common/recurrent artefacts.

- Minor changes:

    - output names of the VCF files.

### `Fixed`

- Contig processing is only available for hg38 reference. ALT and/or HLA contigs can be given external in a file.

    - Automatic generation of HLA/ALT contigs is now possible through tumor BAM instead of fasta.

- Conda links in nf-core modules are fixed.

    - NOTE: Conda environments are not available for the pipeline. Holding conda environment.yml links the same in default creates an error even when enable_conda is false.

- bin/confidenceAnnotation_Indels.py

    - Flag parsing is generic now.
    - Updated to the latest version in https://github.com/DKFZ-ODCF/IndelCallingWorkflow/tree/hg38

- conf/modules.config platypus arguments are fixed (now it is the same as DKFZ/ODCF)

### `Dependencies`

### `Deprecated`
- FREQ-based filtering is removed from bin/confidenceAnnotation_Indels.py.

## v2.0.1 - [09.10.2024]

### `Added`

- nf-prov plugin added

### `Fixed`

- versioning of the pipeline is fixed