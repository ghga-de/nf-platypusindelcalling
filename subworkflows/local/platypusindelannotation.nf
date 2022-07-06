//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF         } from '../../modules/local/annotate_vcf.nf'           addParams( options: params.options )
include { ANNOVAR              } from '../../modules/local/annovar.nf'                addParams( options: params.options )
include { CREATEPIPES          } from '../../modules/local/createpipes.nf'            addParams( options: params.options )
include { CONFIDENCEANNOTATION } from '../../modules/local/confidenceannotation.nf'   addParams( options: params.options )


workflow PLATYPUSINDELANNOTATION {
    take:
    vcf_ch
    sample_ch

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

    versions=Channel.empty()
    // RUN annotate_vcf.pl
    ANNOTATE_VCF (
    vcf_ch, kgenome, dbsnpindel, dbsnpsnv, exac, evs, localcontrolwgs, localcontrolwes, gnomedgenomes, gnomedexomes
    )
    ch_vcf    = ANNOTATE_VCF.out.unziped_vcf
    versions  = versions.mix(ANNOTATE_VCF.out.versions)

    // RUN ANNOVAR and sub-scripts
    ANNOVAR(
    ANNOTATE_VCF.out.forannovar, ch_vcf, annodb
    )
    ch_log   = ANNOVAR.out.log
    versions = versions.mix(ANNOVAR.out.versions)

    // RUN CREATEPIPES (annotate_vcf.pl)
    CREATEPIPES(
    ANNOVAR.out.vcf, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
    )
    versions = versions.mix(CREATEPIPES.out.versions)

    CONFIDENCEANNOTATION(
    CREATEPIPES.out.vcf, sample_ch
    )
    vcf_ch = CONFIDENCEANNOTATION.out.vcf
    conf_vcf_ch   = CONFIDENCEANNOTATION.out.vcf_conf
    versions = versions.mix(CONFIDENCEANNOTATION.out.versions)

emit:
conf_vcf_ch
vcf_ch
versions
}
