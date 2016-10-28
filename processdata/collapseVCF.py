#!/usr/bin/env python

'''
Purpose:       Figure out which SNPs are duplicated from a VCF file and omit the genotype if they are discordant

Usage:         python collapseVCF.py -i infile.vcf -o outfile.vcf
 
Author:        et85, etsai@bwh.harvard.edu

Svn Version:   $Id: collapseVCF.py 2972 2016-01-21 19:32:44Z et85 $

'''

from argparse import ArgumentParser

def inputDuplicateGT(thisGtDict,thisK1,thisK2,thisGT):
    if(thisGT == './.'):
        thisGtDict[thisK1][thisK2][0] += 1
    elif(thisGT == '0/0'):
        thisGtDict[thisK1][thisK2][1] += 1
    elif(thisGT == '0/1'):
        thisGtDict[thisK1][thisK2][2] += 1
    elif(thisGT == '1/1'):
        thisGtDict[thisK1][thisK2][3] += 1
    
    return thisGtDict


### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-i","--infile",
                    type=str,
                    required=False,
                    dest="infile",
                    help="input VCF file")

parser.add_argument("-o","--outfile",
                    type=str,
                    required=False,
                    dest="outfile",
                    help="output VCF file")

args = parser.parse_args()

if(args.infile):
    thisInfile = args.infile

if(args.outfile):
    thisOutfile = args.outfile



thisInFH = open(thisInfile, 'r')
thisOutFH = open(thisOutfile, 'w')

thisDict = {}
thisIdDict = {}

count = 0
prevPos = 0
prevChrPosRefAlt = ''
prevLine = ''
for thisLine in thisInFH.readlines():
    # ignore headers for now
    thisLine = thisLine.rstrip()
    
    if(thisLine.startswith('#CHROM')):
        thisOutFH.write(thisLine + '\n')
        thisHeader = thisLine.split('\t')
        thisSampList = thisHeader[9:]
        nSamp = len(thisHeader[9:])
#         print "nsamp = " + str(nSamp)
    elif(thisLine.startswith("#")):
        thisOutFH.write(thisLine + '\n')
        pass
    else:
        
        count = count + 1
        
        thisRead = thisLine.split('\t')     
        thisChr = thisRead[0]
        thisPos = thisRead[1]
        thisSnpName = thisRead[2]
        thisRef = thisRead[3]
        thisAlt = thisRead[4]
        thisQual = thisRead[5]
        thisFilter = thisRead[6]
        thisInfo = thisRead[7]
        thisFormat = thisRead[8]        
        thisChrPosRefAlt = thisChr + ':' + thisPos + ':' + thisRef + ':' + thisAlt
        
        if(prevChrPosRefAlt != thisChrPosRefAlt): # process the previous variant genotypes (already print prior info)
            if(thisDict.has_key(prevChrPosRefAlt)):
                # there are duplicates
                
                thisSampGtList = []
                for i in range(0,nSamp):        # These keys aren't sorted... fix this
                    # 1 2 == 0, 0 != 0, then 0/0
                    thisSamp = thisSampList[i]
                    if(thisDict[prevChrPosRefAlt][thisSamp][1] != 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
                        and thisDict[prevChrPosRefAlt][thisSamp][3] == 0):
                        thisSampGtList.append('0/0')
                    elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] != 0
                        and thisDict[prevChrPosRefAlt][thisSamp][3] == 0):
                        thisSampGtList.append('0/1')
                    elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
                        and thisDict[prevChrPosRefAlt][thisSamp][3] != 0):
                        thisSampGtList.append('1/1')
                    elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
                        and thisDict[prevChrPosRefAlt][thisSamp][3] == 0 and thisDict[prevChrPosRefAlt][thisSamp][0] != 0):
                        thisSampGtList.append('./.')
#                         print prevChrPosRefAlt + ' ./.'
                    else: # same as above and set to ./. because no consensus
                        thisSampGtList.append('./.')
#                         print prevChrPosRefAlt + ' no consensus'
#                 for thisSamp in (thisDict[prevChrPosRefAlt].keys()):        # These keys aren't sorted... fix this
#                     # 1 2 == 0, 0 != 0, then 0/0
#                     if(thisDict[prevChrPosRefAlt][thisSamp][1] != 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
#                         and thisDict[prevChrPosRefAlt][thisSamp][3] == 0):
#                         thisSampGtList.append('0/0')
#                     elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] != 0
#                         and thisDict[prevChrPosRefAlt][thisSamp][3] == 0):
#                         thisSampGtList.append('0/1')
#                     elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
#                         and thisDict[prevChrPosRefAlt][thisSamp][3] != 0):
#                         thisSampGtList.append('1/1')
#                     elif(thisDict[prevChrPosRefAlt][thisSamp][1] == 0 and thisDict[prevChrPosRefAlt][thisSamp][2] == 0
#                         and thisDict[prevChrPosRefAlt][thisSamp][3] == 0 and thisDict[prevChrPosRefAlt][thisSamp][0] != 0):
#                         thisSampGtList.append('./.')
#                         print prevChrPosRefAlt + ' ./.'
#                     else: # same as above and set to ./. because no consensus
#                         thisSampGtList.append('./.')
#                         print prevChrPosRefAlt + ' no consensus'
                                    
                # print the data set
                
                thisOutFH.write('\t'.join(prevRead[0:2] + [';'.join(thisIdDict[prevChrPosRefAlt])] + prevRead[3:9] + thisSampGtList) + '\n')
#                 print '\t'.join(thisRead[0:2] + [';'.join(thisIdDict[prevChrPosRefAlt])] + thisRead[3:9] + thisSampGtList)
            else:
                ## print previous line -- no duplicates to consider
                if(prevChrPosRefAlt != ''): # skip first line
                    thisOutFH.write(prevLine + '\n')
        else: # else do the following:
#             print "here " + thisChrPosRefAlt
            # save the IDs
            if(prevChrPosRefAlt == thisChrPosRefAlt):
                if(not thisIdDict.has_key(thisChrPosRefAlt)):
                    thisIdDict[thisChrPosRefAlt] = [prevRead[2]]
                
                # Set this read
                thisIdDict[thisChrPosRefAlt].append(thisRead[2])
            
            for i in range(0,nSamp):
                if(not thisDict.has_key(thisChrPosRefAlt)):
                    thisDict[thisChrPosRefAlt] = {}
                
                if(not thisDict[thisChrPosRefAlt].has_key(thisSampList[i])):               
                    thisDict[thisChrPosRefAlt][thisSampList[i]] = [0,0,0,0]

                    # if the previous position was the first one, also need to add it to the count
                    if(prevChrPosRefAlt == thisChrPosRefAlt):
                        thisDict = inputDuplicateGT(thisDict,prevChrPosRefAlt,thisSampList[i],prevRead[i+9])
                
                thisDict = inputDuplicateGT(thisDict,thisChrPosRefAlt,thisSampList[i],thisRead[i+9])
                
#                 if(thisRead[i + 9] == './.'):
#                     thisDict[thisChrPosRefAlt][thisSampList[i]][0] += 1
#                 elif(thisRead[i + 9] == '0/0'):
#                     thisDict[thisChrPosRefAlt][thisSampList[i]][1] += 1
#                 elif(thisRead[i + 9] == '0/1'):
#                     thisDict[thisChrPosRefAlt][thisSampList[i]][2] += 1
#                 elif(thisRead[i + 9] == '1/1'):
#                     thisDict[thisChrPosRefAlt][thisSampList[i]][3] += 1
#             thisDict[thisChrPosRefAlt][thisSampList[i]].append(thisRead[i + 9]) 

        prevPos = thisPos
        prevChrPosRefAlt = thisChrPosRefAlt
        prevRead = thisRead
        prevLine = thisLine

## assuming the last line isn't a duplicate:
if(thisDict.has_key(prevChrPosRefAlt)):
    print "Last line is duplicate"
else:
    thisOutFH.write(prevLine + '\n')


# print thisDict['6:167590750:G:A']['010061321010_R01C01-10007854']



# check dictionary
# for thisK1 in (thisDict.keys()):
#     for thisK2 in (thisDict[thisK1].keys()):
#         # 1 2 == 0, 0 != 0, then 0/0
#         if(thisDict[thisK1][thisK2][1] != 0 and thisDict[thisK1][thisK2][2] == 0
#             and thisDict[thisK1][thisK2][3] == 0):
#             continue
#         elif(thisDict[thisK1][thisK2][1] == 0 and thisDict[thisK1][thisK2][2] != 0
#             and thisDict[thisK1][thisK2][3] == 0):
#             # consensus at 0/1
#             continue
#         elif(thisDict[thisK1][thisK2][1] == 0 and thisDict[thisK1][thisK2][2] == 0
#             and thisDict[thisK1][thisK2][3] != 0):
#             continue
#         elif(thisDict[thisK1][thisK2][1] == 0 and thisDict[thisK1][thisK2][2] == 0
#             and thisDict[thisK1][thisK2][3] == 0 and thisDict[thisK1][thisK2][0] != 0):
#             continue
#         else: # same as above and set to 0/0
#             continue
#         # ./. == ignore for now
#         # 1 2 == 0, 0 != 0, then 0/0
#         # 0 2 == 0, 1 != 0, then 0/1
#         # 0 1 == 0, 2 != 0, then 1/1








