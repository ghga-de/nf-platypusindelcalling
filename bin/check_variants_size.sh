
set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-m max_screenshot]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			somatic_functional_indel_vcf=$2
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

elif [[ $functional_var_count -le $MAX_VARIANT_SCREENSHOTS ]]; then
	printf "Functional variants will be presented\n"
else
    printf "WARNING: No screenshots done, more than $MAX_VARIANT_SCREENSHOTS (cvalue - MAX_VARIANT_SCREENSHOTS) functional variants present in ${somatic_functional_indel_vcf}\n"

fi
