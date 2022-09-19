/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPlatypusindelcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input,
                           params.multiqc_config,
                           params.reference]

def checkPathParamList_annotation = [params.k_genome,
                                     params.dbsnp_indel,
                                     params.dbsnp_snv,
                                     params.exac_file,
                                     params.evs_file,
                                     params.local_control_wgs,
                                     params.local_control_wes,
                                     params.gnomad_genomes,
                                     params.gnomad_exomes,
                                     params.table_folder,
                                     params.annovar_path]

def checkPathParamList_deepanno= [params.repeat_masker,
                                  params.dac_blacklist,
                                  params.duke_excluded,
                                  params.hiseq_depth,
                                  params.self_chain,
                                  params.mapability_file,
                                  params.simple_tandemrepeats,
                                  params.enchancer_file,
                                  params.cpgislands_file,
                                  params.tfbscons_file,
                                  params.encode_dnase_file,
                                  params.mirnas_snornas_file,
                                  params.mirbase_file,
                                  params.cosmic_file,
                                  params.mir_targets_file,
                                  params.cgi_mountains_file,
                                  params.phastconselem_file,
                                  params.encode_tfbs_file]

def checkParamList_runtinda=[params.chrlength_file,
                             params.genemodel_bed,
                             params.exomecapturekit_bed]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.runIndelAnnotation){
    for (param in checkPathParamList_annotation) { if (param) { file(param, checkIfExists: true) } }
    if (params.runIndelDeepAnnotation){
        for (param in checkPathParamList_deepanno) { if (param) { file(param, checkIfExists: true) } }
    }
}

if (params.runTinda){for (param in checkParamList_runtinda) { if (param) { file(param, checkIfExists: true) } }}

//// Check mandatory parameters

// Input samplesheet
if (params.input)       { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// Reference fasta file
if (params.reference)   { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { exit 1, 'Input reference file does not exist' }
// Annotation databases
if (params.k_genome)             { kgenome = Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect() } else { kgenome = Channel.empty() }
if (params.dbsnp_indel)          { dbsnpindel = Channel.fromPath([params.dbsnp_indel, params.dbsnp_indel + '.tbi'], checkIfExists: true).collect() } else { dbsnpindel = Channel.empty() }
if (params.dbsnp_snv)            { dbsnpsnv = Channel.fromPath([params.dbsnp_snv,params.dbsnp_snv +'.tbi' ], checkIfExists: true).collect() } else { dbsnpsnv = Channel.empty() }
if (params.exac_file)            { exac = Channel.fromPath([params.exac_file, params.exac_file + '.tbi'], checkIfExists: true).collect() } else { exac = Channel.empty() }
if (params.evs_file)             { evs = Channel.fromPath([params.evs_file, params.evs_file + '.tbi'], checkIfExists: true).collect() } else { evs = Channel.empty() }
if (params.local_control_wgs)    { localcontrolwgs = Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwgs = Channel.empty() }
if (params.local_control_wes)    { localcontrolwes = Channel.fromPath([params.local_control_wes, params.local_control_wes + '.tbi'], checkIfExists: true).collect() } else { localcontrolwes = Channel.empty() }
if (params.gnomad_genomes)       { gnomadgenomes = Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect() } else { gnomadgenomes = Channel.empty() }
if (params.gnomad_exomes)        { gnomadexomes = Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect() } else { gnomadexomes = Channel.empty() }
// Annovar table folder
if (params.table_folder)         { annodb = Channel.fromPath(params.table_folder, checkIfExists: true ) } else { annodb = Channel.empty() }
// Realiability files
if (params.repeat_masker)        { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.empty() }
if (params.dac_blacklist)        { dacblacklist = Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect() } else { dacblacklist = Channel.empty() }
if (params.duke_excluded)        { dukeexcluded = Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect() } else { dukeexcluded = Channel.empty() }
if (params.hiseq_depth)          { hiseqdepth = Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect() } else { hiseqdepth = Channel.empty() }
if (params.self_chain)           { selfchain = Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect() } else { selfchain = Channel.empty() }
if (params.mapability_file)      { mapability = Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect() } else { mapability = Channel.empty() }
if (params.simple_tandemrepeats) { simpletandemrepeats = Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect() } else { simpletandemrepeats = Channel.empty() }
// Indel Deep Annotation files
if (params.enchancer_file)       { enchangers = Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect() } else { enhangers = Channel.empty() }
if (params.cpgislands_file)      { cpgislands = Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect() } else { cpgislands = Channel.empty() }
if (params.tfbscons_file)        { tfbscons = Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect() } else { tfbscons = Channel.empty() }
if (params.encode_dnase_file)    { encode_dnase = Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect() } else { encode_dnase = Channel.empty() }
if (params.mirnas_snornas_file)  { mirnas_snornas = Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect() } else { mirnas_snornas = Channel.empty() }
if (params.cosmic_file)          { cosmic = Channel.fromPath([params.cosmic_file, params.cosmic_file + '.tbi'], checkIfExists: true).collect() } else { cosmic = Channel.empty() }
if (params.mirbase_file)         { mirbase = Channel.fromPath([params.mirbase_file, params.mirbase_file + '.tbi'], checkIfExists: true).collect() } else { mirbase = Channel.empty() }
if (params.mir_targets_file)     { mir_targets = Channel.fromPath([params.mir_targets_file, params.mir_targets_file + '.tbi'], checkIfExists: true).collect() } else { mir_targets = Channel.empty() }
if (params.cgi_mountains_file)   { cgi_mountains = Channel.fromPath([params.cgi_mountains_file, params.cgi_mountains_file + '.tbi'], checkIfExists: true).collect() } else { cgi_mountains = Channel.empty() }
if (params.phastconselem_file)   { phastconselem = Channel.fromPath([params.phastconselem_file, params.phastconselem_file + '.tbi'], checkIfExists: true).collect() } else { phastconselem = Channel.empty() }
if (params.encode_tfbs_file)     { encode_tfbs = Channel.fromPath([params.encode_tfbs_file, params.encode_tfbs_file + '.tbi'], checkIfExists: true).collect() } else { encode_tfbs = Channel.empty() }
// Tinda files
if (params.chrlength_file)       { chrlength = Channel.fromPath(params.chrlength_file, checkIfExists: true) } else { chrlength = Channel.empty() }
if (params.genemodel_bed)        { genemodel = Channel.fromPath([params.genemodel_bed,params.genemodel_bed +'.tbi'], checkIfExists: true).collect() } else { genemodel = Channel.empty() }

// TODO: Write a pretty log here, write the used parameters
log.info """\
DKFZ-ODCF/IndelCallingWorkflow: A Platypus-based workflow for indel calling
===================================

"""


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include {INPUT_CHECK            } from '../subworkflows/local/input_check'
include {INDEL_CALLING          } from '../subworkflows/local/indel_calling'
include {INDEL_ANNOTATION       } from '../subworkflows/local/indel_annotation'
include {FILTER_VCF             } from '../subworkflows/local/filter_vcf'
include {RUNTINDA               } from '../subworkflows/local/runtinda'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

//  RUN main workflow
workflow PLATYPUSINDELCALLING {

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

//TODO: PREPARE GENOME AND ANNOTATION FILES WITH AWS !!!!

//     SUBWORKFLOW: Read in samplesheet, validate and stage input files

    INPUT_CHECK (
        ch_input
        )

    sample_ch = INPUT_CHECK.out.ch_sample
    sample_ch.view()
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    //
    // SUBWORKFLOW:indelCalling.sh
    //

    INDEL_CALLING(
        sample_ch, ref
    )

    vcf_ch        = INDEL_CALLING.out.vcf_ch
    ch_logs       = ch_logs.mix(INDEL_CALLING.out.ch_platypus_log)
    ch_versions   = ch_versions.mix(INDEL_CALLING.out.platypus_version)

    //
    //SUBWORKFLOW: platypusindelAnnotation.sh
    //
    if (params.runIndelAnnotation) {
        INDEL_ANNOTATION(
        INDEL_CALLING.out.vcf_ch,  kgenome, dbsnpindel, dbsnpsnv, exac, evs, localcontrolwgs, localcontrolwes, gnomadgenomes, gnomadexomes,
        annodb, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats, enchangers,
        cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets, cgi_mountains, phastconselem, encode_tfbs
        )

        ch_versions = ch_versions.mix(INDEL_ANNOTATION.out.versions)
        ch_logs     = ch_logs.mix(INDEL_ANNOTATION.out.logs)
        vcf_ch      = INDEL_ANNOTATION.out.ann_vcf_ch
        conf_vcf    = INDEL_ANNOTATION.out.conf_vcf_ch

        //
        //SUBWORKFLOW: FILTER VCF: filter_vcf.sh
        //
        if (params.runIndelVCFFilter) {
            FILTER_VCF(
            vcf_ch, ref, repeatmasker
            )
            ch_versions = ch_versions.mix(FILTER_VCF.out.versions)
        }
    }
    //
    //SUBWORKFLOW: RUNTINDA: checkSampleSawpTiN.sh
    //

    if (params.runTinda) {
        RUNTINDA(
            INDEL_CALLING.out.vcf_ch, ref, chrlength, genemodel, localcontrolwgs, localcontrolwes, gnomadgenomes, gnomadexomes
            )
        ch_versions = ch_versions.mix(RUNTINDA.out.versions)
        ch_logs     = ch_logs.mix(RUNTINDA.out.logs)
    }

    //
    // MODULE: Pipeline reporting
    //
 //   ch_version_yaml = Channel.empty() 
//    CUSTOM_DUMPSOFTWAREVERSIONS (
 //   ch_versions.unique().collectFile(name: 'collated_versions.yml')
 //   )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPlatypusindelcalling.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
 //  ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_logs.collect{it[1]}
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
