#!/bin/bash

thisOutDir=$1
thisImputeDir=$2
thisSampleListPrefix=$3

# assumption: temp path is created within the anticipated delivery directory
set -e
ls ${thisImputeDir}/*info.gz | xargs -n1 basename | xargs -I {} cp -pr ${thisImputeDir}/{} ${thisOutDir}/${thisSampleListPrefix}_{}
