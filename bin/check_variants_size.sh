
set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-v max_screenshot] [-t tumor] [-c control] [-r ref] [-w window] [-m repeat_masker] [-s iscontrol] [-o screenshots]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			somatic_functional_indel_vcf=$2
			shift # past argument
	    shift # past value
			;;
		-v)
			MAX_VARIANT_SCREENSHOTS=$2
			shift # past argument
	    shift # past value
			;;
		-t)
			FILENAME_TUMOR_BAM=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			FILENAME_CONTROL_BAM=$2
			shift # past argument
	    shift # past value
			;;
		-r)
			REFERENCE_GENOME=$2
			shift # past argument
	    shift # past value
			;;
		-w)
			WINDOW_SIZE=$2
			shift # past argument
	    shift # past value
			;;
		-m)
			REPEAT_MASKER=$2
			shift # past argument
	    shift # past value
			;;
		-s)
			isControlWorkflow=$2
			shift # past argument
	    shift # past value
			;;
		-o)
			combined_screen_shots=$2
			shift # past argument
	    shift # past value
			;;
	esac
done

functional_var_count=`cat ${somatic_functional_indel_vcf} | tail -n +2 | wc -l | cut -f1 -d " "`

if [[ $functional_var_count -eq 0 ]]; then
    printf "WARNING: No functional variants present in ${somatic_functional_indel_vcf}\n"
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

    if [[ "${isControlWorkflow}" == false ]]; then
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

    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${combined_screen_shots} ${sorted}
else
    printf "WARNING: No screenshots done, more than $MAX_VARIANT_SCREENSHOTS (cvalue - MAX_VARIANT_SCREENSHOTS) functional variants present in ${somatic_functional_indel_vcf}\n"

fi
