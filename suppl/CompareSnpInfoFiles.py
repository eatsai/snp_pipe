#!/usr/bin/env python

'''
Purpose:       Take two bpm.csv files, f1 and f2, and output f2 status of f1 probes into a separate file 

Usage:         python compareSnpInfo.py -a file1.snpinfo.txt -b file2.snpinfo.txt
 
Author:        et85, etsai@bwh.harvard.edu
'''

import collections
from argparse import ArgumentParser

####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def crossRef(refDict, queryDict, thisChr, thisPos, thisSnpName):
    if(thisChr in queryDict):
        if(thisPos in queryDict[thisChr]):
            if(thisSnpName in queryDict[thisChr][thisPos]):
                return 1
            else:
                if(sum(len(x) for x in refDict[thisChr][thisPos].itervalues()) > 1):    # if there are more than one probe designed
                    if (thisChr != '0'):
                        return 2
                    else:
                        return 4
                else:
                    return 3
    return 0

####################################################################################
################################    Main Argument   ################################
####################################################################################

### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-a","--file1",
                    type=str,
                    required=False,
                    dest="file1",
                    help="snpinfo.txt input file #1 (file of reference)")

parser.add_argument("-b","--file2",
                    type=str,
                    required=False,
                    dest="file2",
                    help="bpm.csv input file #2 (new file to make comparisons with)")

parser.add_argument("-x","--out1",
                    type=str,
                    required=False,
                    dest="out1",
                    help="output flag file for input file #1 (file of reference)")

parser.add_argument("-y","--out2",
                    type=str,
                    required=False,
                    dest="out2",
                    help="output flag file for input file #2 (new file to make comparisons with)")

args = parser.parse_args()

if(args.file1):
    file1 = args.file1
if(args.file2):
    file2 = args.file2

if(args.out1):
    out1 = args.out1
if(args.out2):
    out2 = args.out2

f1_FH = open(file1, 'r')
f2_FH = open(file2, 'r')

o1_FH = open(out1, 'w')
o2_FH = open(out2, 'w')

file1dict = {}
file2dict = collections.OrderedDict()
# file2dict = {}

# for thisLine in f1_FH.readlines():
#     thisLine = thisLine.rstrip()
#     [thisSnpName, thisChr, thisPos, thisRef, thisAlt, thisRsid, thisEffect] = thisLine.split('\t')
#     
#     file1dict[thisChr][thisPos][thisSnpName] = thisSnpName
#     
# f1_FH.close()


for thisLine in f1_FH.readlines():
    thisLine = thisLine.rstrip()
    [thisSnpName, thisChr, thisPos, thisRef, thisAlt, thisRsid, thisEffect] = thisLine.split('\t')
    
    if(not file1dict.has_key(thisChr)):
        file1dict[thisChr] = collections.OrderedDict()
    if(not file1dict[thisChr].has_key(thisPos)):
        file1dict[thisChr][thisPos] = {}
    
    file1dict[thisChr][thisPos][thisSnpName] = thisSnpName
    
f1_FH.close()



for thisLine in f2_FH.readlines():
    thisLine = thisLine.rstrip()
    [thisSnpName, thisChr, thisPos, thisRef, thisAlt, thisRsid, thisEffect] = thisLine.split('\t')
    
    if(not file2dict.has_key(thisChr)):
        file2dict[thisChr] = collections.OrderedDict()
    if(not file2dict[thisChr].has_key(thisPos)):
        file2dict[thisChr][thisPos] = {}
    
    file2dict[thisChr][thisPos][thisSnpName] = thisSnpName
    
f2_FH.close()


# iterate through the first file dictionary
for thisChr, null in file1dict.items():
    for thisPos, null in file1dict[thisChr].items():
        for thisSnpName, null in file1dict[thisChr][thisPos].items():
            # check each SNP to see if it's in the other dictionary
            thisFlag = crossRef(file1dict, file2dict, thisChr, thisPos, thisSnpName)
            
            o1_FH.write('\t'.join([thisSnpName, thisChr, thisPos, str(thisFlag)]) + '\n')
            
#             print '\t'.join([thisSnpName, thisChr, thisPos, str(thisFlag)])
            
#             print thisChr + ' ' + thisPos + ' ' + thisSnpName


# iterate through the second file dictionary
for thisChr, null in file2dict.items():
    for thisPos, null in file2dict[thisChr].items():
        for thisSnpName, null in file2dict[thisChr][thisPos].items():
            # check each SNP to see if it's in the other dictionary
            thisFlag = crossRef(file2dict, file1dict, thisChr, thisPos, thisSnpName)
            
            o2_FH.write('\t'.join([thisSnpName, thisChr, thisPos, str(thisFlag)]) + '\n')
            
#             print '\t'.join([thisSnpName, thisChr, thisPos, str(thisFlag)])
            
#             print thisChr + ' ' + thisPos + ' ' + thisSnpName


# for thisChr, null in file2dict.items():
#     for thisPos, null in file2dict[thisChr].items():
#         for thisSnpName, null in file2dict[thisChr][thisPos].items():
#             print thisChr + ' ' + thisPos + ' ' + thisSnpName

#         print thisChr + ' ' + thisPos
#         if(thisChr in file2dict):
#             if(thisPos in file2dict[thisChr]):
#                 for thisKey, thisVal in file2dict[thisChr][thisPos]:
#                     print thisChr + ' ' + thisPos + ' ' + thisVal


#                 if(file2dict[thisChr][thisPos].has_key(thisSnpName)):
#                     print thisChr + ' ' + ' ' + thisPos + ' ' + thisSnpName

#         print thisChr + ' ' + thisPos
        









