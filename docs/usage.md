# nf-platypusindelcalling: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Data requirements and Paremeters

Main Swiches:

- --runIndelAnnotation: Run the annotation step or stop the workflow before it.

- --runIndelDeepAnnotation: Run the deep annotation step or stop the workflow before it. Works only of the basic annotation is applied.

- --runIndelVCFFilter: Run the filter step or stop the workflow before it. Works only if the basic annotation is applied.

- --runTinda: Check for sample swaps with TiNDA.

- --skip_multiqc: Skip MultiQC.

Reference Files:

NOTE: this workflow is configured to use both igenomes and refgenie. **genome** can be defined to use either of them. If **genome** parameter is defined then, **fasta**, **fasta_fai**, **chrom_sizes** and **chr_prefix** will be first searched through igenome or refgenie. If not exist those should be manually given.

- --genome: For igenome: GRCh38 or GRCh37, for refgenie: hg38 or hg19

- --fasta: reference fasta file

- --fasta_fai: reference fasta fai file

- --chrom_sizes: Chromosomal size file, should only include primary chromosomes.

- --chr_prefix: chromosomal prefix. tip: left as "" for hg19 reference files, otherwise use "chr"

**Annotation Step:**

NOTE: The reference set bundle which is used in PCAWG study can be found and downloaded [here](https://dcc.icgc.org/api/v1/download?fn=/PCAWG/reference_data/pcawg-dkfz/dkfz-workflow-dependencies_150318_0951.tar.gz). (only in hg19)

If --runIndelAnnotation is true, the following files must be defined (with corresponding indexes):

**1. annotate.vcf Options:**

- --k_genome : 1000k genome indel calls (vcf.gz). Column Name added: 1K_GENOMES
- --dbsnp_indel : dbSNP indel calls (vcf.gz). Column Name added: DBSNP
- --exac_file : ExAC database calls (vcf.gz) (Optional). Column Name added: ExAC
- --evs_file : EVS database calls (vcf.gz) (Optional). Column Name added: EVS
- --local_control_wgs : Extra Local Control WGS file (vcf.gz). Column Name added: LocalControlAF_WGS
- --local_control_wes : Extra Local Control WES file (vcf.gz). Column Name added: LocalControlAF_WES
- --gnomad_genomes : Gnomed Genome sites (vcf.gz). Column Name added: GNOMAD_GENOMES
- --gnomad_exomes : Gnomed Exome sites (vcf.gz). Column Name added: GNOMAD_EXOMES

- --padding : integer
- --minoverlapfraction : float
- --maxborderdist : integer
- --maxmatches : integer

**2. Annovar Options:**

Annovar must be downloaded locally and necessary build version must be generated through:

annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19

- --annovar_path :'/path/to/annovar'

- --buildver : hg19, hg38 or hg18
- --dbtype : wgEncodeGencodeCompV19,wgEncodeGencodeCompV18, wgEncodeGencodeCompV38, wgEncodeGencodeBasicV19, wgEncodeGencodeBasicV18, wgEncodeGencodeBasicV38, refGene, ensGene, knownGene or wgEncodeGencodeBasicV19
- --segdupcol : "SEGDUP_COL"
- --cytobandcol : "CYTOBAND_COL"
- --geneannocols : '"ANNOVAR_FUNCTION,GENE,EXONIC_CLASSIFICATION,ANNOVAR_TRANSCRIPTS"'

**3. Indel Reability Options:**

- --repeat_masker : UCSC Repeat Masker file (bed.gz). Column Name added: REPEAT_MASKER
- --dac_blacklist : UCSC DAC Black List (bed.gz) (Optional). Column Name added: DAC_BLACKLIST
- --duke_excluded : UCSC DUKE Excluded List (bed.gz) (Optional). Column Name added: DUKE_EXCLUDED
- --hiseq_depth : UCSC Hiseq Deep sequencing regions (bed.gz) (Optional). Column Name added: HISEQDEPTH
- --self_chain : UCSC SeLf Chain regions (bed.gz) (Optional). Column Name added: SELFCHAIN
- --mapability_file : UCSC Mappability regions (bed.gz). Column Name added: MAPABILITY
- --simple_tandemrepeats: UCSC Simple tandem repeats (bed.gz). Column Name added: SIMPLE_TANDEMREPEATS

**4. Deep Annotation Options:**

If --runIndelDeepAnnotation is true, at least one of the following files must be defined (with corresponding indexes):

- --enchancer_file : UCSC Enhangers (bed.gz) (Optional). Column Name added: Enhancers
- --cpgislands_file : UCSC CpG islands (bed.gz) (Optional). Column Name added: CpGislands
- --tfbscons_file : UCSC TFBS noncoding sites (bed.gz) (Optional). Column Name added: TFBScons
- --encode_dnase_file : UCSC Encode DNAse cluster (bed.gz) (Optional). Column Name added: ENCODE_DNASE
- --mirna_sncrnas_file : sncRNAs miRBase (bed.gz) (Optional). Column Name added: miRNAs_sncRNAs
- --mirbase_file : miRBase (bed.gz) (Optional). Column Name added: miRBase18
- --cosmic_file : Cosmic coding SNVs (bed.gz) (Optional). Column Name added: COSMIC
- --mir_targets_file : miRNA target sites (bed.gz) (Optional). Column Name added: miRNAtargets
- --cgi_mountains_file : Cgi Mountains (bed.gz) (Optional). Column Name added: CgiMountains
- --phastconselem_file : UCSC Phast Cons Elements (bed.gz) (Optional). Column Name added: phastConsElem20bp
- --encode_tfbs_file : UCSC Encode TFBS (bed.gz) (Optional). Column Name added: ENCODE_TFBS

**Filtration Step:**

If --runVCFFilter is true, the following parameters can be applied:

**1. Filtering Options:**

Filtering is only applied into the samples without control!

**BE CAREFULL** In order to apply below filtartions, annotations to the applied columns must be performed.

- --crit_exac_maxmaf : Max MAF for ExAC. if 0, no filter applys
- --crit_evs_maxmaf : Max MAF for EVC. if 0, no filter applys
- --crit_1kgenomes_maxmaf : Max MAF for 1KGENOMES. if 0, no filter applys
- --crit_localcontrol_maxmaf : Max MAF for LOCAL CONTROL, if 0, no filter applys
- --crit_gnomad_exomes_maxmaf: Max MAF for GNOMAD_GENOMES, if 0, no filter applys
- --crit_gnomad_exomes_maxmaf: Max MAF for GNOMAD_EXOMES, if 0, no filter applys

**2. Indel Extraction Options:**

- --min_confidence_score: integer

**3. Visualize Options:**

- --max_var_screenshots : integer
- --window_size : integer
- --repeat_masker : Repeat Masker file (bed.gz)

**Swap Check Step:**

If --runTinda is true, the following parameters can be applied:

- --genemodel_bed : Genecode Exomes (bed.gz or bed)
- --local_control_wgs : Extra Local Control files (vcf.gz)
- --local_control_wes : Extra Local Control files (vcf.gz)

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

```console
sample,tumor,tumor_index,control,control_index
sample1,sample1_tumor.bam,sample1_tumor.bai,sample1_control.bam,sample1_control.bai
sample2,sample2_tumor.bam,sample2_tumor.bai,sample2_control.bam,sample2_control.bai
sample3,sample3_tumor.bam,sample3_tumor.bai,sample3_control.bam,sample3_control.bai
sample4,sample4_tumor.bam,sample4_tumor.bai,,
sample5,sample5_tumor.bam,sample5_tumor.bai,,
sample6,sample6_tumor.bam,sample6_tumor.bai,,
sample7,sample7_tumor.bam,sample7_tumor.bai,,
```

| Column          | Description                                                                                 |
| --------------- | ------------------------------------------------------------------------------------------- |
| `sample`        | Custom sample name. This entry will be used to name the output files and output directories |
| `tumor`         | Full path to BAM file for tumor samples.                                                    |
| `tumor_index`   | Full path to BAI file for tumor samples.                                                    |
| `control`       | Full path to BAM file for control (optional).                                               |
| `control_index` | Full path to BAI file for control (optional).                                               |

An [example samplesheet](../assets/samplesheet_test.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run main.nf --input samplesheet.csv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

In order to launch the workflow in **DKFZ HPC**, **singularity** must be used with proper config file (text, dkfz_config or dkfz_config_38) which can be found in /conf directory.

```console
module load  nextflow/22.07.1-edge
nextflow run main.nf --input samplesheet_test.csv --outdir results -profile test,singularity
```

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Installation of the pipeline

The pipeline can be pulled from thus github directory and be used!

```console
git clone https://github.com/ghga-de/nf-platypusindelcalling.git
```

Make the bin directory readable:

```console
chmod +x bin/*
```

[Annovar](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) should be downloaded locally if annotation module will be used and should be linked to the appropriare directory in --annovar_path parameter.

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker or Singularity).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`

  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)

- `test`

  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
  - Test only works in dkfz-cluster now

- `dkfz_cluster_38`
  - A profile with a complete configuration for DKFZ cluster
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time.

The computational requirements for the modules can be found and be changed in conf/base.config file

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
