#!/bin/bash

thisSamplistPath=$1
thisOutFilePath=$2
thisBatchVCFPath=$3
thisSnpInfoFilePath=$4

outputVcfFilePath=${thisOutFilePath}.vcf.gz

bcftools view -Oz -S $thisSamplistPath $thisBatchVCFPath -o $outputVcfFilePath
tabix -p vcf $outputVcfFilePath

