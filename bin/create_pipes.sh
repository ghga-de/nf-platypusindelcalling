#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#

#set -x
set -o pipefail


############################################ Create Pipes ####################################################


usage() { echo "Usage: $0 [-i vcf_file] [-id pid_name] [-en enchangers] [-cp cpgislands] [-tf tfbscons] [-ms mirnas_snornas] [-ed encode_dnase] [-mir mirbase] [-c cosmic] [-mt mir_targets] [-cm cgi_mountains] [-p phastconselem] [-et encode_tfbs]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			FILENAME_VCF=$2
			shift # past argument
	    shift # past value
			;;
		-id)
			NAME=$2
			shift # past argument
	    shift # past value
			;;
		-en)
			enchangers=$2
			shift # past argument
	    shift # past value
			;;
		-cp)
			cpgislands=$2
			shift # past argument
	    shift # past value
			;;
		-tf)
			tfbscons=$2
			shift # past argument
	    shift # past value
			;;
		-ms)
			mirnas_snornas=$2
			shift # past argument
	    shift # past value
			;;
		-ed)
			encode_dnase=$2
			shift # past argument
	    shift # past value
			;;
		-mir)
			mirbase=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			cosmic=$2
			shift # past argument
	    shift # past value
			;;
		-mt)
			mir_targets=$2
			shift # past argument
	    shift # past value
			;;
		-cm)
			cgi_mountains=$2
			shift # past argument
	    shift # past value
			;;
		-p)
			phastconselem=$2
			shift # past argument
	    shift # past value
			;;
		-et)
			encode_tfbs=$2
			shift # past argument
	    shift # past value
			;;
	esac
done

if [[ ! -f ${FILENAME_VCF} ]]
then
	echo input file ${FILENAME_VCF} does not exist
	exit 73
fi


output_vcf="${NAME}.deepanno.vcf"

PIPE="zcat < ${FILENAME_VCF}"

if [[ -z  ${enchangers}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${enchangers} --bFileType=bed --columnName='Enhancers'"
fi

if [[ -z  ${cpgislands}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${cpgislands} --bFileType=bed --columnName='CpGislands'"
fi

if [[ -z  ${tfbscons}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${tfbscons} --bFileType=bed --columnName='TFBScons'"
fi

if [[ -z  ${mirnas_snornas}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${mirnas_snornas} --bFileType=bed --columnName='miRNAs_snoRNAs'"
fi

if [[ -z  ${encode_dnase}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${encode_dnase} --bFileType=bed --columnName='ENCODE_DNASE'"
fi

if [[ -z  ${mirbase}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${mirbase} --bFileType=bed --columnName='miRBase18'"
fi

if [[ -z  ${cosmic}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${cosmic} --bFileType=bed --columnName='COSMIC' --bAdditionalColumns=7,8,9 --reportLevel=1"
fi

if [[ -z  ${mir_targets}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${mir_targets} --columnName='miRNAtargets'"
fi

if [[ -z  ${cgi_mountains}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${cgi_mountains} --bFileType=bed --columnName='CgiMountains' --bAdditionalColumns=4"
fi

if [[ -z  ${phastconselem}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${phastconselem} --bFileType=bed --columnName='phastConsElem20bp' --bAdditionalColumns=4"
fi

if [[ -z  ${encode_tfbs}  ]]; then
	PIPE=" ${PIPE} | annotate_vcf.pl -a - -b ${encode_tfbs} --columnName='ENCODE_TFBS'"
fi

eval ${PIPE} > ${output_vcf}

bgzip -c ${output_vcf} > ${output_vcf}.gz

tabix ${output_vcf}.gz