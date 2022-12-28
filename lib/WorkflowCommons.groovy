//
// This file holds several functions common to the multiple workflows in the nf-platypusindelcalling pipeline
//

class WorkflowCommons {

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "============================================================================="
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Function to get column entries from a file
    //
    public static ArrayList getColFromFile(input_file, col=0, uniqify=false, sep='\t') {
        def vals = []
        input_file.eachLine { line ->
            def val = line.split(sep)[col]
            if (uniqify) {
                if (!vals.contains(val)) {
                    vals << val
                }
            } else {
                vals << val
            }
        }
        return vals
    }

    //
    // Function that returns the number of lines in a file
    //
    public static Integer getNumLinesInFile(input_file) {
        def num_lines = 0
        input_file.eachLine { line ->
            num_lines ++
        }
        return num_lines
    }

    //
    // Function to get number of variants reported in BCFTools stats file
    //
    public static Integer getNumVariantsFromBCFToolsStats(bcftools_stats) {
        def num_vars = 0
        bcftools_stats.eachLine { line ->
            def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
            if (matcher) num_vars = matcher[0][1].toInteger()
        }
        return num_vars
    }
}