//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF                  } from '../../modules/local/annotate_vcf.nf'                                 addParams( options: params.options )

workflow PLATYPUSINDELANNOTATION {
    take:
    vcf_ch

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

//exomecapturekitbedfile= file(params.exome_capture_kit_bed_file)
//genemodelbedfile      = file(params.genome_bed_file)

    ANNOTATE_VCF (
    vcf_ch, kgenome, dbsnpindel, dbsnpsnv, exac, evs, localcontrolwgs, localcontrolwes, gnomedgenomes, gnomedexomes
)
    ch_forannovar = ANNOTATE_VCF.out.forannovar
    perl_version = ANNOTATE_VCF.out.versions

emit:
ch_forannovar
perl_version
}
