set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-c isnocontrolcontrol] [-s screenshot] [-r reference] [-n control_bam] [-t tumor_bam] [-w windowsize] [-a repeat_masker] [-m max_screenshot]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			somatic_functional_indel_vcf=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			isNoControlWorkflow=$2
			shift # past argument
	    shift # past value
			;;
		-s)
			combined_screen_shots=$2
			shift # past argument
	    shift # past value
			;;
		-r)
			REFERENCE_GENOME=$2
			shift # past argument
	    shift # past value
			;;
		-n)
			FILENAME_CONTROL_BAM=$2
			shift # past argument
	    shift # past value
			;;
		-t)
			FILENAME_TUMOR_BAM=$2
			shift # past argument
	    shift # past value
			;;
		-w)
			WINDOW_SIZE=$2
			shift # past argument
	    shift # past value
			;;
		-a)
			REPEAT_MASKER=$2
			shift # past argument
	    shift # past value
			;;
		-m)
			MAX_VARIANT_SCREENSHOTS=$2
			shift # past argument
	    shift # past value
			;;
	esac
done

functional_var_count=`cat ${somatic_functional_indel_vcf} | tail -n +2 | wc -l | cut -f1 -d " "`

if [[ $functional_var_count -eq 0 ]]; then
    printf "WARNING: No functional variants present in ${somatic_functional_indel_vcf}\n"
    textPdf.R "WARNING: No functional variants." \
        > "$combined_screen_shots"

elif [[ $functional_var_count -le $MAX_VARIANT_SCREENSHOTS ]]; then

    if [[ "${isControlWorkflow}" == true ]]; then
        visualize.py \
            --vcf=${somatic_functional_indel_vcf} \
            --control=${FILENAME_CONTROL_BAM} \
            --tumor=${FILENAME_TUMOR_BAM} \
            --ref=${REFERENCE_GENOME} \
            --prefix=indel_ \
            --window=${WINDOW_SIZE} \
            --annotations=${REPEAT_MASKER}
    fi

    if [[ "${isNoControlWorkflow}" == true ]]; then
         visualize.py \
            --vcf=${somatic_functional_indel_vcf} \
            --tumor=${FILENAME_TUMOR_BAM} \
            --ref=${REFERENCE_GENOME} \
            --prefix=indel_ \
            --window=${WINDOW_SIZE} \
            --annotations=${REPEAT_MASKER}
    fi
    pngs=(`ls *.pdf`)
    sorted=$(printf "%s\n" ${pngs[@]}|sort -k1,1V)

else
    printf "WARNING: No screenshots done, more than $MAX_VARIANT_SCREENSHOTS (cvalue - MAX_VARIANT_SCREENSHOTS) functional variants present in ${somatic_functional_indel_vcf}\n"
    