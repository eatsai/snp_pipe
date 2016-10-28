#!/bin/bash

thisOutDir=$1
thisOutDirName=$2
thisTempDir=$3
I_readonly=$4
I_compress=$5
I_md5sum=$6
I_debug=$7

if [ $I_readonly = True ]; then
	chmod 440 ${thisOutDir}/*
	chmod 550 ${thisOutDir}
fi

if [ $I_compress = True ]; then
	tar czvf ${thisOutDir}.tgz -C ${thisOutDir}/../ $thisOutDirName
	if [ $I_md5sum = True ]; then
		md5sum ${thisOutDir}.tgz > ${thisOutDir}.tgz.md5
	fi
	chmod 770 ${thisOutDir}
	rm -rf ${thisOutDir}
fi

if [ $I_debug = False ]; then
	rm -r $thisTempDir
fi
