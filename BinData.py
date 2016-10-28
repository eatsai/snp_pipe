#!/usr/bin/env python
import collections
from argparse import ArgumentParser
import math,re,sys

# Purpose:        Take a file and some data in the file. Then, bin a data column
# Input:          1. file name, 2. column with binning data, 3. bin interval
# Example:        
# Author:        et85, etsai@bwh.harvard.edu

parser = ArgumentParser()

parser.add_argument("-f","--file",
                    type=str,
                    dest="infile",
                    required=True,
                    help="path to infile")

parser.add_argument("-o","--outfile",
                    type=str,
                    dest="outfile",
                    required=True,
                    help="path to outfile")

parser.add_argument("-k","--column",
                    type=int,
                    dest="col_num",
                    required=True,
                    help="column number of the data to be binned")

parser.add_argument("-a","--min",
                    type=float,
                    dest="a",
                    default=0,
                    help="min value for binning, each bin will be in the form: a<x<=b")

parser.add_argument("-b","--max",
                    type=float,
                    dest="b",
                    default=1,
                    help="max value for binning, each bin will be in the form: a<x<=b")

parser.add_argument("-i","--int",
                    type=float,
                    dest="bin_int",
                    required=True,
                    help="interval for binning, each bin will be in the form: a<x<=b")

parser.add_argument("-r","--rev_freq",
                    action='store_true',
                    dest="rev_freq",
                    default=False,
                    help="flip frequency values (e.g. new_val = 1-x)")

parser.add_argument("-s","--skip_nheader",
                    type=int,
                    dest="header_rows",
                    default=1,
                    help="number of rows to skip for header [default=1]")

args = parser.parse_args()

def intializeBins(thisMinBinVal,thisMaxBinVal,thisBinInt):
    thisNumBin = int(math.ceil((thisMaxBinVal-thisMinBinVal)/thisBinInt)) + 1 # add one for x=a
    thisBinList = []
    thisBinList.append('['+str(thisMinBinVal)+']')
    for i in range(thisNumBin-1):
        a = thisMinBinVal + thisBinInt * i
        b = thisMinBinVal + thisBinInt * (i+1)
        thisBinList.append('('+str(a)+','+str(b)+']')
    return thisBinList

thisInFH = open(args.infile,"r")
thisOutFH = open(args.outfile,"w")
thisColNum = args.col_num
thisMinBinVal = args.a
thisMaxBinVal = args.b
thisBinInt = args.bin_int
thisNumHeader = args.header_rows

if(args.rev_freq):
    if(thisMinBinVal == 0 and thisMaxBinVal == 1):
        I_RevFreq = args.rev_freq
    else:
        print 'Reverse frequency flag (-r, --rev_freq) is set, but a!=0 or b!=1'
        sys.exit(1)
else:
    I_RevFreq = False

thisBinIntList = intializeBins(thisMinBinVal,thisMaxBinVal,thisBinInt)

thisRow = 0
for thisLine in thisInFH.readlines():
    thisRow += 1
    # squeeze spaces
    thisLine = thisLine.rstrip()
    thisLine = re.sub(' +',' ',thisLine)
    # remove leading/trailing spaces
    thisLine = re.sub('^ ','',thisLine)
    thisLine = re.sub(' $','',thisLine)
    # split by space or tab
    thisRead = re.split(' |\t',thisLine)
    
    # skip headers
    if(thisNumHeader >= thisRow):
        if(I_RevFreq):
            thisOutFH.write('\t'.join(thisRead + ['rev_bin_for_col'+str(thisColNum)]) + '\n')
        else:
            thisOutFH.write('\t'.join(thisRead + ['bin_for_col'+str(thisColNum)]) + '\n')
        continue
    
    if(I_RevFreq):
        thisVal = 1 - float(thisRead[thisColNum - 1])
    else:
        thisVal = float(thisRead[thisColNum - 1])
    
    thisBin = thisBinIntList[int(math.ceil(thisVal/thisBinInt))]
    
    thisOutFH.write('\t'.join(thisRead + [thisBin]) + '\n')
    
#     print thisVal, thisBin
    

thisInFH.close()
thisOutFH.close()
