set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-c isnocontrolcontrol]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			FILENAME_VCF_RAW=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			isNoControlWorkflow=$2
			shift # past argument
	    shift # past value
			;;
	esac
done

# Capture extra VCF columns
totalVCFColumns=11
if [[ ${isNoControlWorkflow} == true ]]; then
  totalVCFColumns=10
fi
extraVCFColumn=$((totalVCFColumns + 1))

corruptLines=`grep -v "^#" ${FILENAME_VCF_RAW} | cut -f ${extraVCFColumn} | sort | uniq | grep -v '^$' | wc -l`
echo "Number of corrupt rows: $corruptLines"

if [[ $corruptLines -gt 0 ]]
then
  (grep "#" ${FILENAME_VCF_RAW} ; grep -v "^#" ${FILENAME_VCF_RAW} | awk -v vcfColumns=$totalVCFColumns '{if(NF == vcfColumns){print $0}}') > ${FILENAME_VCF_RAW}.11
  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 3

  grep -v "^#" ${FILENAME_VCF_RAW} | awk -v vcfColumns=$totalVCFColumns '{if(NF > vcfColumns ){print $0}}' > ${FILENAME_VCF_RAW}.linesCorrupt
  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 5

  if [[ $corruptLines -gt 10 ]]
  then
    printf "Error: Raw VCF should not contain more than 10 columns for the noControl and 11 columns for the control-tumor workflows. The raw file contains more than 10 rows with ${totalVCFColumns} columns.\nIt might have been corrupted or there could be more than one read-groups in a BAM, check ${FILENAME_VCF_RAW}.tmp.platypus.linesCorrupt file\n" && exit 6
  fi
else

  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 2

fi