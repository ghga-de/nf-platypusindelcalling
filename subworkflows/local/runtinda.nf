//
// RUNTINDA: WILL RUN checkSampleSwap_TIN.sh & checkSampleSwap_TIN.pl
//

params.options = [:]

include { SAMPLE_SWAP                                } from '../../modules/local/sample_swap.nf'     addParams( options: params.options )

workflow RUNTINDA {
    take:
    vcf_ch         // channel: [val(meta), vcfgz, vcfgz_tbi]

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }
    if (params.chrlength_file) { chrlength = Channel.fromPath(params.chrlength_file, checkIfExists: true) } else { chrlength = Channel.empty() }
    if (params.genemodel_bed) { genemodel = Channel.fromPath([params.genemodel_bed,params.genemodel_bed +'.tbi'], checkIfExists: true).collect() } else { genemodel = Channel.empty() }
    if (params.local_control_wgs) { localcontrolwgs = Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwgs = Channel.empty() }
    if (params.local_control_wes) { localcontrolwes = Channel.fromPath([params.local_control_wes, params.local_control_wes + '.tbi'], checkIfExists: true).collect() } else { localcontrolwes = Channel.empty() }
    if (params.gnomad_genomes) { gnomadgenomes = Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect() } else { gnomadgenomes = Channel.empty() }
    if (params.gnomad_exomes) { gnomadexomes = Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect() } else { gnomadexomes = Channel.empty() }


    SAMPLE_SWAP(
    vcf_ch, ref, chrlength, genemodel, localcontrolwgs, localcontrolwes, gnomadgenomes, gnomadexomes
    )

    emit:
    vcf_ch
}
