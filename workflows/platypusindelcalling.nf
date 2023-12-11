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
                            params.fasta,
                            params.multiqc_config]

def checkPathParamList_annotation = [params.local_control_wgs,
                                    params.local_control_wes,
                                    params.k_genome,
                                    params.dbsnp_indel,
                                    params.mapability_file,
                                    params.repeat_masker,
                                    params.simple_tandemrepeats]

def checkParamList_runtinda=[params.genemodel_bed,
                            params.local_control_tinda_wgs,
                            params.local_control_tinda_wes,
                            params.gnomad_genomes_tinda,
                            params.gnomad_exomes_tinda]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// To runIndelAnnotation Annovar must be locally installed with its table folder. Other annotation files are required 
// for annotate_vcf.pl 
if (params.runIndelAnnotation){
    for (param in checkPathParamList_annotation) { if (param) { file(param, checkIfExists: true) } }
}

// To runTinda Genemodel bed file, Local control platypus WGS and WES must be provided
if (params.runTinda){
    for (param in checkParamList_runtinda) { if (param) { file(param, checkIfExists: true) } }  
}

// If runIndelDeepAnnotation is true; at least one of the annotation files must be provided
if ((params.runIndelDeepAnnotation) && (!params.enchancer_file && !params.cpgislands_file && !params.tfbscons_file && !params.encode_dnase_file && !params.mirnas_snornas_file && !params.mirna_sncrnas_file && !params.mirbase_file && !params.cosmic_file && !params.mir_targets_file && !params.cgi_mountains_file && !params.phastconselem_file && !params.encode_tfbs_file)) { 
    log.error "Please specify at least one annotation file to perform INDEL Deep Annotation"
    exit 1
}
if (params.annotation_tool.contains("annovar")){
    file(params.annovar_path, checkIfExists: true)
}

//// Check mandatory parameters

// Reference genome
ref            = Channel.fromPath([params.fasta,params.fasta_fai], checkIfExists: true).collect()
chr_prefix     = Channel.value(params.chr_prefix)
chrlength      = params.chrom_sizes ? Channel.fromPath(params.chrom_sizes, checkIfExists: true) : Channel.empty()   

// Input samplesheet
if (params.input)                { ch_input  = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Annotation databases
kgenome             =  params.k_genome            ? Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
dbsnpindel          =  params.dbsnp_indel         ? Channel.fromPath([params.dbsnp_indel, params.dbsnp_indel + '.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
exac                = params.exac_file            ? Channel.fromPath([params.exac_file, params.exac_file + '.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
evs                 = params.evs_file             ? Channel.fromPath([params.evs_file, params.evs_file + '.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
localcontrolwgs     = params.local_control_wgs    ? Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
localcontrolwes     = params.local_control_wes    ? Channel.fromPath([params.local_control_wes,params.local_control_wes + '.tbi' ], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]}     
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
gnomadgenomes       = params.gnomad_genomes       ? Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])
gnomadexomes        = params.gnomad_exomes        ? Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect().map{[[id: "annotate_vcf_ref"], it[0],it[1]]}     
                                                  : Channel.of([[id: "annotate_vcf_ref"],[],[]])

annotate_vcf_ref    = kgenome.join(dbsnpindel)
                            .join(exac)
                            .join(evs)
                            .join(localcontrolwgs)
                            .join(localcontrolwes)
                            .join(gnomadgenomes)
                            .join(gnomadexomes)
                                                  
// Annovar table folder
annodb              = params.annovar_path         ? Channel.fromPath(params.annovar_path + '/humandb/') 
                                                  : Channel.empty()
// VEP cache
vep_cache_db        = params.vep_cache          ? Channel.fromPath(params.vep_cache).collect()         : []

// Realiability files
repeatmasker        = params.repeat_masker        ? Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
dacblacklist        = params.dac_blacklist        ? Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
dukeexcluded        = params.duke_excluded        ? Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
hiseqdepth          = params.hiseq_depth          ? Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
selfchain           = params.self_chain           ? Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
mapability          = params.mapability_file      ? Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])
simpletandemrepeats = params.simple_tandemrepeats ? Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect().map{[[id: "realibility_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "realibility_ref"],[],[]])

realibility_ref     = repeatmasker.join(dacblacklist) 
                                .join(dukeexcluded) 
                                .join(hiseqdepth)
                                .join(selfchain)
                                .join(mapability)
                                .join(simpletandemrepeats)    
                                                                       
// Indel Deep Annotation files
enchangers          = params.enchancer_file       ? Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
cpgislands          = params.cpgislands_file      ? Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
tfbscons            = params.tfbscons_file        ? Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
encode_dnase        = params.encode_dnase_file    ? Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
mirnas_snornas      = params.mirnas_snornas_file  ? Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]}
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
mirnas_sncrnas      = params.mirna_sncrnas_file   ? Channel.fromPath([params.mirna_sncrnas_file, params.mirna_sncrnas_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
cosmic              = params.cosmic_file          ? Channel.fromPath([params.cosmic_file, params.cosmic_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
mirbase             = params.mirbase_file         ? Channel.fromPath([params.mirbase_file, params.mirbase_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
mir_targets         = params.mir_targets_file     ? Channel.fromPath([params.mir_targets_file, params.mir_targets_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
cgi_mountains       = params.cgi_mountains_file   ? Channel.fromPath([params.cgi_mountains_file, params.cgi_mountains_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
phastconselem       = params.phastconselem_file   ? Channel.fromPath([params.phastconselem_file, params.phastconselem_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]} 
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
encode_tfbs         = params.encode_tfbs_file     ? Channel.fromPath([params.encode_tfbs_file, params.encode_tfbs_file + '.tbi'], checkIfExists: true).collect().map{[[id: "deepanno_ref"], it[0],it[1]]}
                                                  : Channel.of([[id: "deepanno_ref"],[],[]])
deepanno_ref        = enchangers.join(cpgislands)
                                .join(tfbscons)
                                .join(encode_dnase)
                                .join(mirnas_snornas)
                                .join(cosmic)
                                .join(mirbase)
                                .join(mir_targets)
                                .join(cgi_mountains)
                                .join(phastconselem)
                                .join(encode_tfbs)
                                .join(mirnas_sncrnas)
                                                 
// Tinda files
genemodel            = params.genemodel_bed           ? Channel.fromPath(params.genemodel_bed, checkIfExists: true).collect() 
                                                      : Channel.empty()
localcontroltindawgs = params.local_control_tinda_wgs ? Channel.fromPath([params.local_control_tinda_wgs,params.local_control_tinda_wgs + '.tbi' ], checkIfExists: true).collect().map{[[id: "tinda_ref"], it[0],it[1]]} 
                                                      : Channel.of([[id: "tinda_ref"],[],[]])
localcontroltindawes = params.local_control_tinda_wes ? Channel.fromPath([params.local_control_tinda_wes, params.local_control_tinda_wes + '.tbi'], checkIfExists: true).collect().map{[[id: "tinda_ref"], it[0],it[1]]} 
                                                      : Channel.of([[id: "tinda_ref"],[],[]])
gnomadgenomes_tinda  = params.gnomad_genomes_tinda    ? Channel.fromPath([params.gnomad_genomes_tinda, params.gnomad_genomes_tinda + '.tbi'], checkIfExists: true).collect().map{[[id: "tinda_ref"], it[0],it[1]]}  
                                                      : Channel.of([[id: "tinda_ref"],[],[]])
gnomadexomes_tinda   = params.gnomad_exomes_tinda     ? Channel.fromPath([params.gnomad_exomes_tinda, params.gnomad_exomes_tinda + '.tbi'], checkIfExists: true).collect().map{[[id: "tinda_ref"], it[0],it[1]]} 
                                                      : Channel.of([[id: "tinda_ref"],[],[]])
tinda_ref           = localcontroltindawgs.join(localcontroltindawes)
                                            .join(gnomadgenomes_tinda)
                                            .join(gnomadexomes_tinda)

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

//
// MODULE: Local Modules
//

include { GREP_SAMPLENAME       } from '../modules/local/grep_samplename.nf'
include { SAMPLE_SWAP           } from '../modules/local/sample_swap.nf'
include { GETCHROMSIZES         } from '../modules/local/getchromsizes.nf'

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

//  RUN main workflow
workflow PLATYPUSINDELCALLING {

    // To gather all versions.yml for MultiQC
    ch_versions = Channel.empty()
    //To gather all logs for MultiQC
    ch_logs = Channel.empty()

    //    
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //    
    INPUT_CHECK (
        ch_input
        )
    sample_ch = INPUT_CHECK.out.ch_sample
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Extract sample name from BAM
    //
    GREP_SAMPLENAME(
        sample_ch
        )
    ch_versions = ch_versions.mix(GREP_SAMPLENAME.out.versions)

    if ( !params.chrom_sizes) {
        //
        // MODULE: Prepare chromosome size file if not provided
        //
        GETCHROMSIZES(
            ref
            )
        ch_versions = ch_versions.mix(GETCHROMSIZES.out.versions)
        chrlength   = GETCHROMSIZES.out.sizes
    }

    //
    // SUBWORKFLOW:indelCalling.sh
    //   
    INDEL_CALLING(
        sample_ch, 
        ref
    )
    ch_logs         = ch_logs.mix(INDEL_CALLING.out.log_ch)
    ch_versions     = ch_versions.mix(INDEL_CALLING.out.versions)

    // Prepare an input channel of vcf with sample names 
    INDEL_CALLING.out.vcf_ch
        .join(GREP_SAMPLENAME.out.samplenames)
        .set{ch_vcf}

    //
    //SUBWORKFLOW: platypusindelAnnotation.sh
    //
    // annotation has two part, first annotation for annovar, second is deep annotation includes for various genomic regions 
    if (params.runIndelAnnotation) {
        INDEL_ANNOTATION(
            ch_vcf, 
            annotate_vcf_ref, 
            realibility_ref, 
            deepanno_ref, 
            chr_prefix,
            ref,
            annodb,
            vep_cache_db 
        )
        ch_versions = ch_versions.mix(INDEL_ANNOTATION.out.versions)

        //
        //SUBWORKFLOW: FILTER VCF: filter_vcf.sh
        //
        // filtering is only apply into the samples without control, 
        // indel extraction and visualization will apply both cases.
        // runIndelVCFFilter is only applicable if basic annotation done.  
        
        if (params.runIndelVCFFilter) {
            FILTER_VCF(
                INDEL_ANNOTATION.out.ann_vcf_ch, 
                ref, 
                repeatmasker
            )
            ch_versions = ch_versions.mix(FILTER_VCF.out.versions)
        }
        else{
            println "Skipping indel vcf filtering"
        }
    }
    else{
        println "Skipping indel annotation"
    }
    //
    // MODULE: RUNTINDA
    //
    // Checks sample swap in platypus output vcf
    // checkSampleSawpTiN.sh
    if (params.runTinda) {
        SAMPLE_SWAP(
            ch_vcf, 
            ref, 
            chrlength, 
            genemodel, 
            tinda_ref, 
            chr_prefix
            )
        ch_versions = ch_versions.mix(SAMPLE_SWAP.out.versions)    
        ch_logs     = ch_versions.mix(SAMPLE_SWAP.out.log) 
    }
    else{
        println "Skipping sample swap check"
    }

    // Info required for completion email and summary
    def multiqc_report = []
    if (!params.skip_multiqc){
        //
        // MODULE: Pipeline reporting
        //
        ch_version_yaml = Channel.empty() 
        CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'))

        //
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowPlatypusindelcalling.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_logs.collect{it[1]}.ifEmpty({[]}))

        MULTIQC (
            ch_multiqc_files.collect()
            )
        
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }
    else{
        println "Skipping MultiQC"
    }
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
