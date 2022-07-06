//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )


workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK (samplesheet)
        .csv
        .splitCsv ( header:true, sep:',' )
        .map{ create_bam_channel(it) }
        .set {ch_sample}

    emit:
    ch_sample // channel: [ sample, tumor,tumor.bai, control, control.bai, iscontrol ]
    versions = SAMPLESHEET_CHECK.out.versions
}

// Function to get list of [ sample, [ tumor, control ] ]
def create_bam_channel(LinkedHashMap row) {

    // add path(s) of the fastq file(s) to the meta map
    def bam_meta = []
        if (!file(row.tumor).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Tumor file does not exist!\n${row.tumor}"
        }

        if (row.iscontrol) {
            if (!file(row.control).exists()) {
                if (row.control == 'dummy') {bam_meta = [  row.sample, file(row.tumor), file(row.tumor + '.bai'),[],[], row.iscontrol  ]}
                else {exit 1, "ERROR: Please check input samplesheet -> Control file does not exist!\n${row.control}"}
            }
            bam_meta = [ row.sample, file(row.tumor), file(row.tumor + '.bai'), file(row.control), file(row.control + '.bai'),row.iscontrol ]

        } else {
            bam_meta = [  row.sample, file(row.tumor), file(row.tumor + '.bai'),[],[], row.iscontrol  ]
        }
    return bam_meta
}


