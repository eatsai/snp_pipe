#!/bin/bash

thisInfilePath=$1
thisOutFilePath=$2
thisSnpInfoFilePath=$3

mergeVcfFilePath="${thisOutFilePath}.vcf.gz"

bcftools merge -l $thisInfilePath -O z -o $mergeVcfFilePath
tabix -p vcf ${mergeVcfFilePath}

if [ -f $thisSnpInfoFilePath ]; then
	cp -p $thisSnpInfoFilePath ${thisOutFilePath}.snpinfo.txt
else
	echo Directory does not exist: $thisSnpInfoFilePath
	exit 1
fi
