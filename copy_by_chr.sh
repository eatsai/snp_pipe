#!/bin/bash

thisFile=$1
thisImputeDir=$2
thisOutPath=$3
thisVcfName=chr${LSB_JOBINDEX}.dose.vcf

# assumption: temp path is created within the anticipated delivery directory
cp -p ${thisImputeDir}/${thisVcfName}.gz ${thisOutPath}/${thisFile}.${thisVcfName}.gz

tabix -p vcf ${thisOutPath}/${thisFile}.${thisVcfName}.gz
