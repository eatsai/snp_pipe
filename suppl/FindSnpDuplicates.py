#!/usr/bin/env python

## 2015-01-08 Temporarily shelved. This code may not be necessary.

'''
Purpose:       Figure out which SNPs are duplicated from a given map file

Usage:         python FindSnpDuplicates.py -m infile.map > outfile.txt
 
Author:        et85, etsai@bwh.harvard.edu
'''

from argparse import ArgumentParser

### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-m","--map",
                    type=str,
                    required=False,
                    dest="map",
                    help="input map file")

parser.add_argument("-a","--alleleinfo",
                    type=str,
                    required=False,
                    dest="alleleinfo",
                    help="input allele info file")

args = parser.parse_args()

if(args.map):
    thisMapFile = args.map

if(args.alleleinfo):
    thisAlleleInfoFile = args.alleleinfo

thisSnpDict = {}

###  This is for map file
thisMapFileFH = open(thisMapFile, "r")
for thisLine in thisMapFileFH.readlines():
    thisLine = thisLine.rstrip()
    (thisChr, thisSnpName, thisCM, thisPos) = thisLine.split("\t")
    if(thisChr != '0'):
        thisChrPos = thisChr + ':' + thisPos
        if(thisSnpDict.has_key(thisChrPos)):
            thisSnpDict[thisChrPos].append(thisSnpName)
        else:
            thisSnpDict[thisChrPos] = [thisSnpName]

# thisAlleleInfoFileFH = open(thisAlleleInfoFile, "r")
# for thisLine in thisAlleleInfoFileFH.readlines():
#     thisLine = thisLine.rstrip()
#     (thisSnpName, thisChr, thisPos, thisFwdRef, thisFwdAlt, thisStrand) = thisLine.split("\t")
#     if(thisChr != '0'):
#         thisKeyName = thisChr + ':' + thisPos + ':' + thisFwdRef + ':' + thisFwdAlt
#         if(thisSnpDict.has_key(thisKeyName)):
#             thisSnpDict[thisKeyName].append(thisSnpName)
#         else:
#             thisSnpDict[thisKeyName] = [thisSnpName]

count = 0
thisSum = 0
for thisKey in thisSnpDict.keys():
    if (len(thisSnpDict[thisKey]) > 1):
        thisSum = thisSum + len(thisSnpDict[thisKey])
#         count = count + 1
#         print thisKey,    
#         print thisSnpDict[thisKey]
#     
#     if (count > 10):
#         break

print thisSum
    
    





