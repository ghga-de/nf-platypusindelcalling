//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF as ANNOTATE_VCF1    } from '../../modules/local/annotate_vcf.nf'                       addParams( options: params.options )
include { ANNOVAR as ANNOVAR1              } from '../../modules/local/annovar.nf'                            addParams( options: params.options )
include { CREATEPIPES as CREATEPIPES1      } from '../../modules/local/createpipes.nf'                        addParams( options: params.options )
include { ANNOTATE_VCF as ANNOTATE_VCF2    } from '../../modules/local/annotate_vcf.nf'                       addParams( options: params.options )
include { ANNOVAR as ANNOVAR2              } from '../../modules/local/annovar.nf'                            addParams( options: params.options )
include { CREATEPIPES as CREATEPIPES2      } from '../../modules/local/createpipes.nf'                        addParams( options: params.options )
include { CONFIDENCEANNOTATION_WITHCONTROL } from '../../modules/local/confidenceannotation_withcontrol.nf'   addParams( options: params.options )
include { CONFIDENCEANNOTATION_NOCONTROL   } from '../../modules/local/confidenceannotation_nocontrol.nf'     addParams( options: params.options )


workflow PLATYPUSINDELANNOTATION {
    take:
    vcf_ch_withcontrol
    vcf_ch_nocontrol
    sample_withcontrol
    sample_nocontrol

    main:

    if (params.k_genome) { kgenome = Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect() } else { kgenome = Channel.empty() }
    if (params.dbsnp_indel) { dbsnpindel = Channel.fromPath([params.dbsnp_indel, params.dbsnp_indel + '.tbi'], checkIfExists: true).collect() } else { dbsnpindel = Channel.empty() }
    if (params.dbsnp_snv) { dbsnpsnv = Channel.fromPath([params.dbsnp_snv,params.dbsnp_snv +'.tbi' ], checkIfExists: true).collect() } else { dbsnpsnv = Channel.empty() }
    if (params.exac_file) { exac = Channel.fromPath([params.exac_file, params.exac_file + '.tbi'], checkIfExists: true).collect() } else { exac = Channel.empty() }
    if (params.evs_file) { evs = Channel.fromPath([params.evs_file, params.evs_file + '.tbi'], checkIfExists: true).collect() } else { evs = Channel.empty() }
    if (params.local_control_wgs) { localcontrolwgs = Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwgs = Channel.empty() }
    if (params.local_control_wes) { localcontrolwes = Channel.fromPath([params.local_control_wes, params.local_control_wes + '.tbi'], checkIfExists: true).collect() } else { localcontrolwes = Channel.empty() }
    if (params.gnomed_genomes) { gnomedgenomes = Channel.fromPath([params.gnomed_genomes, params.gnomed_genomes + '.tbi'], checkIfExists: true).collect() } else { gnomedgenomes = Channel.empty() }
    if (params.gnomed_exomes) { gnomedexomes = Channel.fromPath([params.gnomed_exomes, params.gnomed_exomes + '.tbi'], checkIfExists: true).collect() } else { gnomedexomes = Channel.empty() }

    if (params.repeat_masker) { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.empty() }
    if (params.dac_blacklist) { dacblacklist = Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect() } else { dacblacklist = Channel.empty() }
    if (params.duke_excluded) { dukeexcluded = Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect() } else { dukeexcluded = Channel.empty() }
    if (params.hiseq_depth) { hiseqdepth = Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect() } else { hiseqdepth = Channel.empty() }
    if (params.self_chain) { selfchain = Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect() } else { selfchain = Channel.empty() }
    if (params.mapability_file) { mapability = Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect() } else { mapability = Channel.empty() }
    if (params.simple_tandemrepeats) { simpletandemrepeats = Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect() } else { simpletandemrepeats = Channel.empty() }
    if (params.table_folder) { annodb = Channel.fromPath(params.table_folder, checkIfExists: true) } else { annodb = Channel.empty() }


    withcontrol_vcf_ch = Channel.empty()
    nocontrol_vcf_ch = Channel.empty()
    perl_version = Channel.empty()


///// PIPELINE RUN WITH CONTROL SAMPLES /////

    if (vcf_ch_withcontrol){
    // RUN annotate_vcf.pl
        ANNOTATE_VCF1 (
        vcf_ch_withcontrol, kgenome, dbsnpindel, dbsnpsnv, exac, evs, localcontrolwgs, localcontrolwes, gnomedgenomes, gnomedexomes
        )
        ch_forannovar = ANNOTATE_VCF1.out.forannovar
        ch_vcf        = ANNOTATE_VCF1.out.unziped_vcf
        perl_version  = ANNOTATE_VCF1.out.versions

    // RUN ANNOVAR and sub-scripts
        ANNOVAR1(
        ch_forannovar, ch_vcf, annodb
        )
        ch_log = ANNOVAR1.out.log
    //    annovar_version=ANNOVAR.out.versions

    // RUN CREATEPIPES (annotate_vcf.pl)
        CREATEPIPES1(
        ANNOVAR1.out.vcf, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
        )

        CONFIDENCEANNOTATION_WITHCONTROL(
        CREATEPIPES1.out.vcf, sample_withcontrol
        )
        withcontrol_vcf_ch =CONFIDENCEANNOTATION_WITHCONTROL.out.vcf
    }

///// PIPELINE RUN WITH NO CONTROL SAMPLES /////
    if (vcf_ch_nocontrol){
        // RUN annotate_vcf.pl
        ANNOTATE_VCF2 (
        vcf_ch_nocontrol, kgenome, dbsnpindel, dbsnpsnv, exac, evs, localcontrolwgs, localcontrolwes, gnomedgenomes, gnomedexomes
        )
        ch_forannovar = ANNOTATE_VCF2.out.forannovar
        ch_vcf        = ANNOTATE_VCF2.out.unziped_vcf
        perl_version  = ANNOTATE_VCF2.out.versions

        // RUN ANNOVAR and sub-scripts
        ANNOVAR2(
        ch_forannovar, ch_vcf, annodb
        )
        ch_log = ANNOVAR2.out.log
        //    annovar_version=ANNOVAR.out.versions

        // RUN CREATEPIPES (annotate_vcf.pl)
        CREATEPIPES2(
        ANNOVAR2.out.vcf, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
        )
        CONFIDENCEANNOTATION_NOCONTROL(
        CREATEPIPES2.out.vcf, sample_nocontrol
        )

        nocontrol_vcf_ch=CONFIDENCEANNOTATION_NOCONTROL.out.vcf
    }


emit:
withcontrol_vcf_ch
nocontrol_vcf_ch
perl_version
}
