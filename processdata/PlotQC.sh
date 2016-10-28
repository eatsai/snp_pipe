#!/bin/bash

# Purpose: This script will trigger a generation of the daily QC plots

export PYTHONPATH=/usr/lib/python2.7/site-packages:$PYTHONPATH

thisDate=`date --date="1 days ago" +%Y%m%d`

/usr/local/bin/python /PHShome/et85/git/biobank/MakeQCPlots.py -B /PHShome/et85/biobank/ref/productfiles/Multi-EthnicGlobal_A1.csv -O /mnt/et85/phs-pcpgm-research/ResCore_ProductionData/GTP_Production/QCPlots/$thisDate -D /mnt/illumina/illumina_lims/call/$thisDate

####### set the cron job to run every day at 4am:
# $ ssh ngs_pcpgm.dipr.partners.org
# $ crontab -e
# * 4 * * * /PHShome/et85/bin/PlotSnpQC.sh >/dev/null 2>&1

## Please note that >/dev/null 2>&1 prevents the cron job from sending emails.