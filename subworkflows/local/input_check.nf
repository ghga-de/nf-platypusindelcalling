//
// Check input samplesheet and get read channels
//

params.options = [:]

//include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .map{it ->
            def sample = it.sample
            def tumor = file(it.tumor)
            def tumor_bai = file(it.tumor + '.bai')
            def control = file(it.control)
            def control_bai = file(it.control + '.bai')
            return [sample, tumor, tumor_bai, control, control_bai]
            }
        .set { ch_sample }
    emit:
    ch_sample // channel: [ meta, tumor, control ]
}

