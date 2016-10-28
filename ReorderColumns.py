#!/usr/bin/env python
import collections
from argparse import ArgumentParser
import math,re,sys

# Purpose:       Reorder columns
# Input:
# Author:        et85, etsai@bwh.harvard.edu
parser = ArgumentParser()


parser.add_argument("-i","--infile",
                    type=str,
                    dest="infile",
                    required=True,
                    help="path to infile")

parser.add_argument("-o","--outfile",
                    type=str,
                    dest="outfile",
                    required=True,
                    help="path to file2")

parser.add_argument("-k","--column",
                    type=int,
                    nargs='+',
                    dest="col_list",
                    required=True,
                    help="columns to display in desired order")

args = parser.parse_args()

thisInFH = open(args.infile, "r")
thisOutFH = open(args.outfile, "w")
thisColList = args.col_list

# print thisColList

for thisLine in thisInFH.readlines():
    # squeeze spaces
    thisLine = thisLine.rstrip()
    thisLine = re.sub(' +',' ',thisLine)
    # remove leading/trailing spaces
    thisLine = re.sub('^ ','',thisLine)
    thisLine = re.sub(' $','',thisLine)
    # split by space or tab
    thisRead = re.split(' |\t',thisLine)
    
    thisOutFH.write('\t'.join(thisRead[col_num-1] for col_num in thisColList) + '\n')
#     print '\t'.join(thisRead[col_num-1] for col_num in thisColList)
    
    
    
thisInFH.close()
thisOutFH.close()











