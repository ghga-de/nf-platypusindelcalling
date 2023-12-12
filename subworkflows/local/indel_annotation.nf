//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF           } from '../../modules/local/annotate_vcf.nf'                     addParams( options: params.options )
include { ANNOVAR                } from '../../modules/local/annovar.nf'                          addParams( options: params.options )
include { INDEL_RELIABILITY_PIPE } from '../../modules/local/indel_reliability_pipe.nf'           addParams( options: params.options )
include { CONFIDENCE_ANNOTATION  } from '../../modules/local/confidence_annotation.nf'            addParams( options: params.options )
include { ANNOTATION_PIPES       } from '../../modules/local/annotation_pipes.nf'                 addParams( options: params.options )
include { ENSEMBLVEP_VEP         } from '../../modules/nf-core/modules/ensemblvep/vep/main'       addParams( options: params.options )
include { ENSEMBLVEP_DOWNLOAD    } from '../../modules/nf-core/modules/ensemblvep/download/main'  addParams( options: params.options )


workflow INDEL_ANNOTATION {
    take:
    vcf_ch           // channel: [val(meta), vcf.gz, vcf.gz.tbi ,val(tumorname), val(controlname) ]
    kgenome          // channel: [file,index]
    dbsnpindel       // channel: [file,index]
    exac             // channel: [file,index]
    evs              // channel: [file,index]
    localcontrolwgs  // channel: [file,index] 
    localcontrolwes  // channel: [file,index]
    gnomadgenomes    // channel: [file,index]
    gnomadexomes     // channel: [file,index]
    repeatmasker     // channel: [file,index]
    dacblacklist     // channel: [file,index]
    dukeexcluded     // channel: [file,index]
    hiseqdepth       // channel: [file,index]
    selfchain        // channel: [file,index]
    mapability       // channel: [file,index]
    simpletandemrepeats  // channel: [file,index]
    enchangers       // channel: [file,index]
    cpgislands       // channel: [file,index]
    tfbscons         // channel: [file,index]
    encode_dnase     // channel: [file,index] 
    mirnas_snornas   // channel: [file,index]
    cosmic           // channel: [file,index]
    mirbase          // channel: [file,index]
    mir_targets      // channel: [file,index]
    cgi_mountains    // channel: [file,index]
    phastconselem    // channel: [file,index]
    encode_tfbs      // channel: [file,index]
    mirnas_sncrnas   // channel: [file,index]
    chr_prefix       // val channel: [prefix]
    ref              // channel [fasta,fai] 
    annodb           // channel: [table_annovar_dir]
    vep_cache        // channel: [vep_cache_dir]

    main:

    versions=Channel.empty()
    logs=Channel.empty() 
    //
    // MODULE: ANNOTATE_VCF
    //
    // RUN annotate_vcf.pl: Uses various databases (all mandatory) to annotate variants
    input_ch = vcf_ch.map{ it -> tuple( it[0], it[1], it[2])}
    ANNOTATE_VCF (
        input_ch, 
        kgenome,dbsnpindel,exac,evs,localcontrolwgs,localcontrolwes,gnomadgenomes,gnomadexomes, 
        chr_prefix
    )
    versions  = versions.mix(ANNOTATE_VCF.out.versions)

    ANNOTATE_VCF.out.unziped_vcf
        .join(ANNOTATE_VCF.out.forannovar)
        .set{ input_ch}

    if (params.annotation_tool.contains("annovar")){
        //
        // MODULE: ANNOVAR
        //
        // RUN annovar, processAnnovarOutput.pl and newCols2vcf.pl: annovar annotates and classifies the variants, 
        // perl scripts re-creates vcfs.
        ANNOVAR(
            input_ch, 
            annodb, 
            chr_prefix
        )
        logs          = logs.mix(ANNOVAR.out.log)
        versions      = versions.mix(ANNOVAR.out.versions)
        annotated_vcf = ANNOVAR.out.vcf
        }
    else{
        
        if(params.download_cache){
            ENSEMBLVEP_DOWNLOAD(
                input_ch.map{ it -> tuple( it[0], it[1])}
                )
            versions  = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
            vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
        }
        
        ENSEMBLVEP_VEP(
            ANNOTATE_VCF.out.unziped_vcf,
            vep_cache,
            ref
        )
        versions      = versions.mix(ENSEMBLVEP_VEP.out.versions)
        annotated_vcf = ENSEMBLVEP_VEP.out.vcf
    }

    //
    // MODULE: INDEL_RELIABILITY_PIPE
    //
    // RUN annotate_vcf.pl : BED files are used to annotate variants
    INDEL_RELIABILITY_PIPE(
        annotated_vcf, 
        repeatmasker,dacblacklist,dukeexcluded,hiseqdepth,selfchain,mapability,simpletandemrepeats
    )
    versions = versions.mix(INDEL_RELIABILITY_PIPE.out.versions)

    //
    // MODULE: CONFIDENCE_ANNOTATION
    //
    // RUN: confidenceAnnotation_Indels.py : Confidence annotation will be added to the variants
    if (params.fasta.contains("38")){
        ref_type = "hg38"   
    }
    else{
        ref_type = "hg37"
    }
    input_ch = vcf_ch.join(INDEL_RELIABILITY_PIPE.out.vcf)
    input_ch = input_ch.map{ it -> tuple( it[0], it[3], it[4], it[5], it[6])}
    CONFIDENCE_ANNOTATION(
        input_ch, 
        ref_type
    )
    ann_vcf_ch  = CONFIDENCE_ANNOTATION.out.vcf_ann
    versions    = versions.mix(CONFIDENCE_ANNOTATION.out.versions)

    // RUN annotate_vcf.pl : Uses optional databases to annotate variants, only given databases will be used. 
    if (params.runIndelDeepAnnotation)
    {
        //
        // MODULE: ANNOTATION_PIPES
        //
        ANNOTATION_PIPES (
            ann_vcf_ch, 
            enchangers, cpgislands,tfbscons,encode_dnase,mirnas_snornas,cosmic,mirbase,mir_targets,cgi_mountains,phastconselem,encode_tfbs,mirnas_sncrnas
        )
        ann_vcf_ch  = ANNOTATION_PIPES.out.vcf 
        versions    = versions.mix(ANNOTATION_PIPES.out.versions)
    }
    else{
        println "Skipping deep annotation since runIndelDeepAnnotation is set to false"
    }
emit:
logs
ann_vcf_ch
versions
}
