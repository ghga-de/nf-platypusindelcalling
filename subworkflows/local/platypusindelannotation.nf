//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF                  } from '../../modules/local/annotate_vcf.nf'                                 addParams( options: params.options )

workflow PLATYPUSINDELANNOTATION {
    take:
    vcf_ch

    main:

    if (params.k_genome) { kgenome = Channel.fromPath(params.k_genome, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.dbsnp_indel) { dbsnpindel = Channel.fromPath(params.dbsnp_indel, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.dbsnp_snv) { dbsnpsnv = Channel.fromPath(params.dbsnp_snv, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.exac_file) { exac = Channel.fromPath(params.exac_file, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.evs_file) { evs = Channel.fromPath(params.evs_file, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.local_control_wgs) { localcontrolwgs = Channel.fromPath(params.local_control_wgs, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.local_control_wes) { localcontrolwes = Channel.fromPath(params.local_control_wes, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.gnomed_genomes) { gnomedgenomes = Channel.fromPath(params.gnomed_genomes, checkIfExists: true) } else { ref = Channel.empty() }
    if (params.gnomed_exomes) { gnomedexomes = Channel.fromPath(params.gnomed_exomes, checkIfExists: true) } else { ref = Channel.empty() }

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
