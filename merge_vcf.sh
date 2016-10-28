#!/bin/bash

thisSampMergeListPath=$1
thisOutFilePath=$2
thisSnpInfoFilePath=$3

outputVcfFilePath=${thisOutFilePath}.vcf.gz

bcftools merge -l $thisSampMergeListPath -O z -o $outputVcfFilePath
tabix -p vcf $outputVcfFilePath

if [ -f $thisSnpInfoFilePath ]; then
	cp -p $thisSnpInfoFilePath ${thisOutFilePath}.snpinfo.txt
else
	echo Directory does not exist: $thisSnpInfoFilePath
	exit 1
fi

