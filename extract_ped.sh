#!/bin/bash

thisPedPathList=$1
thisOutFilePath=$2
thisMapFilePath=$3
thisSnpinfoFilePath=$4

outputVcfFilePath=${thisOutFilePath}.ped

touch $outputVcfFilePath
for thisPedPath in `cat $thisPedPathList`; do
	cat $thisPedPath >> $outputVcfFilePath
done

# copy map files
if [ -f $thisMapFilePath ]; then
	cp -p $thisMapFilePath ${thisOutFilePath}.map
else
	echo Directory does not exist: $thisMapFilePath
	exit 1
fi

# copy snpinfo files
if [ -f $thisSnpinfoFilePath ]; then
	cp -p $thisSnpinfoFilePath ${thisOutFilePath}.snpinfo.txt
else
	echo Directory does not exist: $thisSnpinfoFilePath
	exit 1
fi
