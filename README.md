

[![GitHub Actions CI Status](https://github.com/kubranarci/nf-core-platypusindelcalling/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/platypusindelcalling/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/kubranarci/nf-core-platypusindelcalling/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/platypusindelcalling/actions?query=workflow%3A%22nf-core+linting%22)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)



## Introduction

**DKFZ-ODCF/nf-platypusindelcalling**:A Platypus-based insertion/deletion-detection workflow with extensive quality control additions.  The workflow is based on DKFZ - ODCF OTP Indel Calling Pipeline.

For now, this workflow is only optimal to work in ODCF Cluster. The config file (conf/dkfz_cluster.config) can be used as an example. Running Annotation, DeepAnnotation, Filter and Tinda steps are optinal and can be turned off using [runIndelAnnotation, runIndelDeepAnnotation, runIndelVCFFilter, runTinda] parameters sequentialy.  

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

The Indel workflow was in the pan-cancer analysis of whole genomes (PCAWG) and can be cited with the following publication:

Pan-cancer analysis of whole genomes.
The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium.
Nature volume 578, pages 82–93 (2020).
DOI 10.1038/s41586-020-1969-6


## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Platypus ([`Platypus`](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data))
   : Platypus tool is used to call variants using local realignmnets and local assemblies. It can detect SNPs, MNPs, short indels, replacements, deletions up to several kb. It can be both used with WGS and WES. The tool has been thoroughly tested on data mapped with Stampy and BWA.
2. Basic Annotations: In-house scripts to annotate with several databases like gnomAD, dbSNP, and ExAC.
3. ANNOVAR ([`Annovar`](https://annovar.openbioinformatics.org/en/latest/))
   : annotate_variation.pl is used to annotate variants. The tool makes classifications for intergenic, intogenic, nonsynoymous SNP, frameshift deletion or large-scale duplication regions.
4. Reliability and confidation annotations: It is an optional step (--runIndelAnnotation) for Mapability, hiseq, selfchain and repeat regions checks for reliability and confidence of those scores.
5. INDEL Deep Annotation: It is an optional step (runIndelDeepAnnotation) for number of extra indel annotations like enhancer, cosmic, mirBASE, encode databases.
6. Filtering: It is an optional step (runIndelVCFFilter)for only applies for the tumor samples with no-control.
7. Checks Sample Swap: Canopy Based Clustering and Bias Filter

## Data requirements and Paremeters

Reference Files:

- --reference: A reference directory .fa or .fasta (including index files in the same directory) should be defined.

**Platypus Options**

Platypus parameters must be given under platypus_params, fow whole list of options check assets/platypus_options.txt

- --platypus_params    : " --genIndels=1 --genSNPs=1 --verbosity=1 --bufferSize=100000 --maxReads=5000000 --minFlank=0"


**Annotation Step:** 

If --runIndelAnnotation is true, the following files must be defined (with corresponding indexes):

**1. annotate.vcf Options:**

- --k_genome           : 1000k genome indel calls (vcf.gz)
- --dbsnp_indel        : dbSNP indel calls (vcf.gz)
- --dbsnp_snv          : dbSNP SV calls (vcf.gz)
- --exac_file          : ExAC database calls (vcf.gz)
- --evs_file           : EVS database calls (vcf.gz)
- --local_control_wgs  : Extra Local Control files (vcf.gz)
- --local_control_wes  : Extra Local Control files (vcf.gz) 
- --gnomad_genomes     : Gnomed Genome sites (vcf.gz)
- --gnomad_exomes      : Gnomed Exome sites (vcf.gz)

- --padding             : integer
- --minoverlapfraction : float
- --maxborderdist       : integer
- --maxmatches          : integer

**2. Annovar Options:**

Annovar must be downloaded locally and necessary build version must be generated through:

annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19

- --table_folder        :'/path/to/annovar/humandb/'
- --annovar_path        :'/path/to/annovar'

- --buildver            : hg19, hg38 or hg18
- --dbtype              : wgEncodeGencodeCompV19,wgEncodeGencodeCompV18, wgEncodeGencodeCompV38, wgEncodeGencodeBasicV19, wgEncodeGencodeBasicV18, wgEncodeGencodeBasicV38, refGene, ensGene, knownGene, wgEncodeGencodeBasicV19
- --segdupcol           : "SEGDUP_COL"
- --cytobandcol         : "CYTOBAND_COL"
- --geneannocols        : '"ANNOVAR_FUNCTION,GENE,EXONIC_CLASSIFICATION,ANNOVAR_TRANSCRIPTS"'

**3. Indel Reability Options:**

- --repeat_masker       : Repeat Masker file (bed.gz)
- --dac_blacklist       : Black List (bed.gz)
- --duke_excluded       : Excluded List (bed.gz)
- --hiseq_depth         : Hiseq Deep sequencing regions (bed.gz)
- --self_chain          : SeLf Chain regions (bed.gz)
- --mapability_file     : Mappability regions (bed.gz)
- --simple_tandemrepeats: Simple tandem repeats (bed.gz)

**4. Confidence Annotation Options:**

Confidence annotation parameters must be given under confidence_opts_indel, fow whole list of options check assets/confidenceAnnotation_Indel_options.txt

- --confidence_opts_indel: " --hetcontr=-4.60517 --tumaltgen=0" 

**5. Deep Annotation Options:**

If --runIndelDeepAnnotation is true, the following files must be defined (with corresponding indexes):

-  --enchancer_file      : Enhangers (bed.gz)
-  --cpgislands_file     : CpG islands (bed.gz)
-  --tfbscons_file       : TFBS noncoding sites (bed.gz)
-  --encode_dnase_file   : Encode DNAse cluster (bed.gz)
-  --mirnas_snornas_file : snoRNAs miRBase  (bed.gz)
-  --mirbase_file        : miRBase (bed.gz)
-  --cosmic_file         : Cosmic coding SNVs (bed.gz)
-  --mir_targets_file    : miRNA target sites (bed.gz)
-  --cgi_mountains_file  : Cgi Mountains (bed.gz)
-  --phastconselem_file  : Phast Cons Elements (bed.gz)
-  --encode_tfbs_file    :  Encode TFBS (bed.gz)

**Filtration Step:** 

If --runVCFFilter is true, the following parameters can be applied:

**1. Filtering Options:**

Filtering is only applied into the samples without control! Filtering options must be inserted into filter_values parameter. 

Available filtering columns if the databases are provided: 
add " EVS MAF VALUE+" to filter_values parameter for EVS db (1)

add " ExAC AF VALUE+" to filter_values parameter for ExAC db (0.01)

add " GNOMAD_EXOMES AF VALUE+" to filter_values parameterfor gnomAD exomes db (0.001)

add " GNOMAD_GENOMES AF VALUE+" to filter_values parameter for gnomAD genomes db (0.001)

add " DBSNP CLN,COMMON nonexist,exist" to filter_values parameter to filtering clinical variants

add " 1K_GENOMES EUR_AF VALUE+" to filter_values parameter for 1k genomes db  for EUR (0.01)

add " LocalControlAF_WGS AF VALUE+ LocalControlAF_WES AF VALUE+" to filter_values parameter for Local control (0.01)

add " LocalControlAF_WGS . VALUE+" to filter_values parameter to filter recurrance (7)

- --filter_values       : " GNOMAD_EXOMES AF 0.001+ GNOMAD_GENOMES AF 0.001+"

**2. Indel Extraction Options:**

- --min_confidence_score: integer

**3. Visualize Options:**

- --max_var_screenshots : integer
- --window_size         : integer
- --repeat_masker       : Repeat Masker file (bed.gz)

**Swap Check Step:** 

If --runTinda is true, the following parameters can be applied:

- -- chrlength_file     : 'assets/chrlengths.plain.sorted.txt'
- --genemodel_bed       :  Genecode Exomes (bed.gz)
- --exomecapturekit_bed : Exome capture kit UTR regions (bed.gz)
- --local_control_wgs   : Extra Local Control files (vcf.gz)
- --local_control_wes   : Extra Local Control files (vcf.gz) 
- --gnomad_genomes      : Gnomed Genome sites (vcf.gz)
- --gnomad_exomes       : Gnomed Exome sites (vcf.gz)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/))

3. Download the pipeline and test it on a minimal dataset with a single command:

   git clone https://github.com/kubranarci/nf-platypusindelcalling.git

  before run do this to bin directory, make it runnable!:

  ```console
  chmod +x bin/*
  ```

   ```console
   nextflow run nf-core/platypusindelcalling -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   
4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run nf-core/platypusindelcalling --input samplesheet.csv --outdir <OUTDIR> --profile <docker/singularity/institute>
   ```
## Samplesheet columns

**sample**: The sample name will be tagged to the job

**tumor**: The path to the tumor file

**control**: The path to the control file, if there is no control will be kept blank.

## Documentation

**TODO**
The nf-core/platypusindelcalling pipeline comes with documentation about the pipeline [usage](https://nf-co.re/platypusindelcalling/usage), [parameters](https://nf-co.re/platypusindelcalling/parameters) and [output](https://nf-co.re/platypusindelcalling/output).


## Credits

nf-core/platypusindelcalling was originally written by Kuebra Narci kuebra.narci@dkfz-heidelberg.de.

The pipeline is originally written in workflow management language Roddy. [Inspired github page] (https://github.com/DKFZ-ODCF/IndelCallingWorkflow)

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#platypusindelcalling` channel](https://nfcore.slack.com/channels/platypusindelcalling) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/platypusindelcalling for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:
