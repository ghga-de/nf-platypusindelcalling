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

def checkParamList_runtinda=[params.genemodel_bed,
                            params.local_control_platypus_wgs,
                            params.local_control_platypus_wes]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// To runIndelAnnotation Annovar must be locally installed with its table folder. Other annotation files are required 
// for annotate_vcf.pl 
if (params.runIndelAnnotation){
    for (param in checkPathParamList_annotation) { if (param) { file(param, checkIfExists: true) } }
}
// To runTinda Genemodel bed file, Local control platypus WGS and WES must be provided
if (params.runTinda){for (param in checkParamList_runtinda) { if (param) { file(param, checkIfExists: true) } }}

// If runIndelDeepAnnotation is true; at least one of the annotation files must be provided
if ((params.runIndelDeepAnnotation) && (!params.enchancer_file && !params.cpgislands_file && !params.tfbscons_file && !params.encode_dnase_file && !params.mirnas_snornas_file && !params.mirbase_file && !params.cosmic_file && !params.mir_targets_file && !params.cgi_mountains_file && !params.phastconselem_file && !params.encode_tfbs_file)) { 
    log.error "Please specify at least one annotation file to perform INDEL Deep Annotation"
    exit 1
}

//// Check mandatory parameters

// Input samplesheet
if (params.input)                { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
if (params.annovar_path)         { annodb = Channel.fromPath(params.annovar_path + '/humandb/', checkIfExists: true ) } else { annodb = Channel.empty() }
// Realiability files
if (params.repeat_masker)        { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.empty() }
if (params.dac_blacklist)        { dacblacklist = Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect() } else { dacblacklist = Channel.empty() }
if (params.duke_excluded)        { dukeexcluded = Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect() } else { dukeexcluded = Channel.empty() }
if (params.hiseq_depth)          { hiseqdepth = Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect() } else { hiseqdepth = Channel.empty() }
if (params.self_chain)           { selfchain = Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect() } else { selfchain = Channel.empty() }
if (params.mapability_file)      { mapability = Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect() } else { mapability = Channel.empty() }
if (params.simple_tandemrepeats) { simpletandemrepeats = Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect() } else { simpletandemrepeats = Channel.empty() }
// Indel Deep Annotation files
if (params.enchancer_file)       { enchangers = Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect() } else { enchangers = Channel.of([],[]) }
if (params.cpgislands_file)      { cpgislands = Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect() } else { cpgislands = Channel.of([],[]) }
if (params.tfbscons_file)        { tfbscons = Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect() } else { tfbscons = Channel.of([],[]) }
if (params.encode_dnase_file)    { encode_dnase = Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect() } else { encode_dnase = Channel.of([],[]) }
if (params.mirnas_snornas_file)  { mirnas_snornas = Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect() } else { mirnas_snornas = Channel.of([],[]) }
if (params.cosmic_file)          { cosmic = Channel.fromPath([params.cosmic_file, params.cosmic_file + '.tbi'], checkIfExists: true).collect() } else { cosmic = Channel.of([],[]) }
if (params.mirbase_file)         { mirbase = Channel.fromPath([params.mirbase_file, params.mirbase_file + '.tbi'], checkIfExists: true).collect() } else { mirbase = Channel.of([],[]) }
if (params.mir_targets_file)     { mir_targets = Channel.fromPath([params.mir_targets_file, params.mir_targets_file + '.tbi'], checkIfExists: true).collect() } else { mir_targets = Channel.of([],[]) }
if (params.cgi_mountains_file)   { cgi_mountains = Channel.fromPath([params.cgi_mountains_file, params.cgi_mountains_file + '.tbi'], checkIfExists: true).collect() } else { cgi_mountains = Channel.of([],[]) }
if (params.phastconselem_file)   { phastconselem = Channel.fromPath([params.phastconselem_file, params.phastconselem_file + '.tbi'], checkIfExists: true).collect() } else { phastconselem = Channel.of([],[]) }
if (params.encode_tfbs_file)     { encode_tfbs = Channel.fromPath([params.encode_tfbs_file, params.encode_tfbs_file + '.tbi'], checkIfExists: true).collect() } else { encode_tfbs = Channel.of([],[]) }
// Tinda files
if (params.genemodel_bed)        { genemodel = Channel.fromPath([params.genemodel_bed,params.genemodel_bed +'.tbi'], checkIfExists: true).collect() } else { genemodel = Channel.empty() }
if (params.local_control_platypus_wgs)    { localcontrolplatypuswgs = Channel.fromPath([params.local_control_platypus_wgs,params.local_control_platypus_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolplatypuswgs = Channel.empty() }
if (params.local_control_platypus_wes)    { localcontrolplatypuswes = Channel.fromPath([params.local_control_platypus_wes, params.local_control_platypus_wes + '.tbi'], checkIfExists: true).collect() } else { localcontrolplatypuswes = Channel.empty() }

// Save AWS IGenomes file containing annotation version
//def anno_readme = params.genomes[params.genome]?.readme
//if (anno_readme && file(anno_readme).exists()) {
//    file("${params.outdir}/genome/").mkdirs()
//    file(anno_readme).copyTo("${params.outdir}/genome/")
//}

// Set up reference depending on the genome choice
// NOTE: link will be defined by aoutomatic reference generation when the pipeline ready!
if (params.ref_type)
    {
    if (params.ref_type == 'hg37')
        { 
        def fa_file  = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef_Phix/hs37d5_PhiX.fa"
        ref          = Channel.fromPath([fa_file,fa_file +'.fai'], checkIfExists: true).collect() 
        def chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hs37d5.fa.chrLenOnlyACGT_realChromosomes.tab'
        chrlength    = Channel.fromPath(chr_file, checkIfExists: true)
        chr_prefix   = Channel.value("")
        }
    if (params.ref_type == 'hg19') 
        { 
        def fa_file = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_1-22_X_Y_M.fa"
        ref = Channel.fromPath([fa_file,fa_file +'.fai'], checkIfExists: true).collect()
        def chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_1-22_X_Y_M.fa.chrLenOnlyACGT.tab'
        chrlength = Channel.fromPath(chr_file, checkIfExists: true)
        chr_prefix   = Channel.value("chr")  
        }
    }
else
{
    if (params.reference)      { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { exit 1, 'Input reference file does not exist' }
    if (params.chrlength_file) { chrlength = Channel.fromPath(params.chrlength_file, checkIfExists: true) } else { exit 1, 'Chromosome length file does not exist'  }
    if (params.chr_prefix)     {chr_prefix= Channel.of(params.chr_prefix)} else {chr_prefix= Channel.value("")}
}

// TODO: Write a pretty log here, write the used parameters
log.info """\
DKFZ-ODCF/IndelCallingWorkflow: A Platypus-based workflow for indel calling
===================================
    input                 : ${params.input}
    outdir                : ${params.outdir}
    ref_type              : ${params.ref_type}

    Post Processes
    runIndelAnnotation    : ${params.runIndelAnnotation}
    runIndelDeepAnnotation: ${params.runIndelDeepAnnotation}
    runIndelVCFFilter     : ${params.runIndelVCFFilter}
    runTinda              : ${params.runTinda}
    skip_multiqc          : ${params.skip_multiqc}

    Annotate VCF
    padding               : ${params.padding}
    minoverlapfraction    : ${params.minoverlapfraction} 
    maxborderdist         : ${params.maxborderdist} 
    maxmatches            : ${params.maxmatches}

    Annovar
    annovar_path          : ${params.annovar_path}
    buildver              : ${params.buildver}
    dbtype                : ${params.dbtype}
    cytobandcol           : ${params.cytobandcol}
    segdupcol             : ${params.segdupcol}
    geneannocols          : ${params.geneannocols} 

    Filtering
    min_confidence_score  : ${params.min_confidence_score}
    filter_exac           : ${params.filter_exac}
    filter_1kgenomes      : ${params.filter_1kgenomes}
    filter_recurrance     : ${params.filter_recurrance}
    filter_localcontrol   : ${params.filter_localcontrol}
    filter_non_clinic     : ${params.filter_non_clinic}
    crit_exac_maxmaf      : ${params.crit_exac_maxmaf}
    crit_evs_maxmaf       : ${params.filter_exac}
    crit_1kgenomes_maxmaf : ${params.crit_1kgenomes_maxmaf}
    crit_recurrance       : ${params.crit_recurrance}
    crit_localcontrol_maxmaf: ${params.crit_localcontrol_maxmaf}

    Visualization
    max_var_screenshots   : ${params.max_var_screenshots}
    window_size           : ${params.window_size}

    Tinda
    seqtype               : ${params.seqtype}
    right_border          : ${params.right_border}
    bottom_border         : ${params.bottom_border}
    maf_threshold         : ${params.maf_threshold}
"""
.stripIndent()

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

//include {SET_CHR            } from '../modules/local/set_chr.nf'
include { GREP_SAMPLENAME   } from '../modules/local/grep_samplename.nf'
include { SAMPLE_SWAP      } from '../modules/local/sample_swap.nf'

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
    //
    // SUBWORKFLOW:indelCalling.sh
    //
    
    INDEL_CALLING(
        sample_ch, ref
    )

    ch_logs         = ch_logs.mix(INDEL_CALLING.out.ch_platypus_log)
    ch_versions     = ch_versions.mix(INDEL_CALLING.out.platypus_version)

    // Prepare an input channel of vcf with sample names 
    name_ch=GREP_SAMPLENAME.out.samplenames
    vcf_ch=INDEL_CALLING.out.vcf_ch
    ch_vcf=vcf_ch.join(name_ch)

    //
    //SUBWORKFLOW: platypusindelAnnotation.sh
    //
    // annotation has two part, first annotation for annovar, second is deep annotation includes for various genomic regions 
    if (params.runIndelAnnotation) {
        INDEL_ANNOTATION(
        ch_vcf, kgenome, dbsnpindel, exac, evs, localcontrolwgs,
        localcontrolwes, gnomadgenomes, gnomadexomes, annodb, repeatmasker, dacblacklist,
        dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats, enchangers,
        cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets,
        cgi_mountains, phastconselem, encode_tfbs, chr_prefix 
        )

        ch_versions = ch_versions.mix(INDEL_ANNOTATION.out.versions)

        //
        //SUBWORKFLOW: FILTER VCF: filter_vcf.sh
        //
        // filtering is only apply into the samples without control, 
        // indel extraction and visualization will apply both cases.
        // runIndelVCFFilter is only applicable if annotation done.  
        
        if (params.runIndelVCFFilter) {
            FILTER_VCF(
            INDEL_ANNOTATION.out.ann_vcf_ch, ref, repeatmasker
            )
            ch_versions = ch_versions.mix(FILTER_VCF.out.versions)
        }
    }
    //
    //SUBWORKFLOW: RUNTINDA: checkSampleSawpTiN.sh
    //
    // Checks sample swap in platypus output vcf
    if (params.runTinda) {

        SAMPLE_SWAP(
            ch_vcf, ref, chrlength, genemodel, localcontrolplatypuswgs, 
            localcontrolplatypuswes, gnomadgenomes, gnomadexomes, chr_prefix
            )
        ch_versions = ch_versions.mix(SAMPLE_SWAP.out.versions)    
        ch_logs = ch_versions.mix(SAMPLE_SWAP.out.log) 
    }

    // Info required for completion email and summary
    def multiqc_report = []
    ch_versions.view()
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
