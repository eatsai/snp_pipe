#!/bin/bash

thisFile=$1
thisImputeDir=$2
thisOutPath=$3
thisTempPath=$4
thisVcfName=chr${LSB_JOBINDEX}.dose.vcf

# assumption: temp path is created within the anticipated delivery directory
vcftools --keep ${thisTempPath}/${thisFile}.txt --gzvcf ${thisImputeDir}/${thisVcfName}.gz --recode --stdout | bgzip > ${thisOutPath}/${thisFile}.${thisVcfName}.gz

tabix -p vcf ${thisOutPath}/${thisFile}.${thisVcfName}.gz
