/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPlatypusindelcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.reference]
def checkPathParamList_annotation = [params.k_genome, params.dbsnp_indel, params.dbsnp_snv,params.exac_file, params.evs_file,
                                     params.local_control_wgs, params.local_control_wes, params.gnomed_genomes, params.gnomed_exomes,
                                     params.table_folder, params.annovar_path]

def checkPathParamList_deepanno= [params.repeat_masker, params.dac_blacklist, params.duke_excluded, params.hiseq_depth,
                                  params.self_chain, params.mapability_file, params.simple_tandemrepeats]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.runIndelAnnotation){
    for (param in checkPathParamList_annotation) { if (param) { file(param, checkIfExists: true) } }
    if (params.runIndelDeepAnnotation){
        for (param in checkPathParamList_deepanno) { if (param) { file(param, checkIfExists: true) } }
    }
}

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


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

//
// MODULE: Local modules
//
include { EXTRACT_SAMPLE_NAME } from '../modules/local/extract_sample_name.nf'     addParams( options: params.options )

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

// MODULE: Prepare sample names
    EXTRACT_SAMPLE_NAME(
        sample_ch
        )

    //
    // SUBWORKFLOW:indelCalling.sh
    //

    INDEL_CALLING(sample_ch)

    platypus_vcf_ch = INDEL_CALLING.out.vcf_ch
    ch_logs         = ch_logs.mix(INDEL_CALLING.out.ch_platypus_log)
    ch_versions     = ch_versions.mix(INDEL_CALLING.out.platypus_version)

    name_tumor=EXTRACT_SAMPLE_NAME.out.tumor_name
    name_control=EXTRACT_SAMPLE_NAME.out.control_name
    //
    //SUBWORKFLOW: platypusindelAnnotation.sh
    //
    if (params.runIndelAnnotation) {
        INDEL_ANNOTATION(
        platypus_vcf_ch, name_tumor, name_control
        )
        ch_versions = ch_versions.mix(INDEL_ANNOTATION.out.versions)
        vcf_ch      = INDEL_ANNOTATION.out.vcf_ch
        conf_vcf    = INDEL_ANNOTATION.out.conf_vcf_ch

        //
        //SUBWORKFLOW: FILTER VCF: filter_vcf.sh
        //
        if (params.runIndelVCFFilter) {
            FILTER_VCF(
            vcf_ch
            )
            ch_versions = ch_versions.mix(FILTER_VCF.out.versions)
        }
    }
    //
    //SUBWORKFLOW: RUNTINDA: checkSampleSawpTiN.sh
    //

    if (params.runTinda) {
        RUNTINDA(
            platypus_vcf_ch
            )
    }

    //
    // MODULE: Pipeline reporting
    //
   // CUSTOM_DUMPSOFTWAREVERSIONS (
  //  ch_versions.unique().collectFile(name: 'collated_versions.yml')
  //  )


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPlatypusindelcalling.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
   // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

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
