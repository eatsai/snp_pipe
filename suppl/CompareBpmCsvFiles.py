#!/usr/bin/env python

'''
Purpose:       Take two bpm.csv files, f1 and f2, and compare which one has 

Usage:         python compareBpmCsvFiles.py -a file1.bpm.csv -b file2.bpm.csv
 
Author:        et85, etsai@bwh.harvard.edu
'''

from argparse import ArgumentParser

####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def compareData(thisDict, thisIlmnID, thisSnpName, thisIlmnStrand, thisTopGenSeq):
    [dictSnpName, dictIlmnStrand, dictTopGenSeq] = thisDict[thisIlmnID]
    
#     print [thisSnpName, dictSnpName, thisIlmnStrand, dictIlmnStrand, thisTopGenSeq, dictTopGenSeq]
    # check to make sure everything aligns
    thisFlag = 0
    if(dictSnpName != thisSnpName):
        thisFlag += 1
    if(dictIlmnStrand != thisIlmnStrand):
        thisFlag += 1
    if(dictTopGenSeq != thisTopGenSeq):
        thisFlag += 1
    
    return thisFlag
    
    

####################################################################################
################################    Main Argument   ################################
####################################################################################

### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-a","--file1",
                    type=str,
                    required=False,
                    dest="file1",
                    help="bpm.csv input file #1 (file of reference)")

parser.add_argument("-b","--file2",
                    type=str,
                    required=False,
                    dest="file2",
                    help="bpm.csv input file #2 (new file to make comparisons with)")

parser.add_argument("-o","--outfile",
                    type=str,
                    required=False,
                    dest="outfile",
                    help="output file with illumina flags for SNP array in file #2 with indicators of probes found in file #1")

parser.add_argument("-f","--flagname",
                    type=str,
                    required=False,
                    dest="flagname",
                    help="name of the flag file")

args = parser.parse_args()

if(args.file1):
    file1 = args.file1

if(args.file2):
    file2 = args.file2

if(args.outfile):
    outfile = args.outfile
    I_out = 1
else:
    I_out = 0

if(args.flagname):
    flagname = args.flagname

thisF1Dict = {}
thisDict2 = {} # may not need this

count = 0
countOmitted = 0
countDifferent = 0
countOK = 0
countDiffProbe = 0

with open(file1, 'r') as thisFile1FH:
    for _ in xrange(8):  #skips the header range
        next(thisFile1FH)
    for thisLine in thisFile1FH:
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split(',')

        if thisLine.find('[Controls]') != -1:
            break
        
        thisIlmnID = thisRead[0]
        thisSnpName = thisRead[1]
        thisIlmnStrand = thisRead[2]
        thisTopGenSeq = thisRead[17]
        
        thisF1Dict[thisIlmnID] = [thisSnpName, thisIlmnStrand, thisTopGenSeq]
        thisDict2[thisSnpName] = [thisIlmnID, thisIlmnStrand, thisTopGenSeq]
    
#         print [thisIlmnID, thisSnpName, thisIlmnStrand, thisTopGenSeq]        
#         print thisLine
#         count += 1
#         if(count > 10):
#             break

thisFile1FH.close()

# print output
outFH = open(outfile, 'w')
outFH.write('Name\t' + flagname + '\n')

with open(file2, 'r') as thisFile2FH:
    for _ in xrange(8):  #skips the header range
        next(thisFile2FH)
    for thisLine in thisFile2FH:
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split(',')

        if thisLine.find('[Controls]') != -1:
            break
        
        thisIlmnID = thisRead[0]
        thisSnpName = thisRead[1]
        thisIlmnStrand = thisRead[2]
        thisTopGenSeq = thisRead[17]
        I_probe = 0
        
        if(thisF1Dict.has_key(thisIlmnID)):
            thisFlag = compareData(thisF1Dict, thisIlmnID, thisSnpName, thisIlmnStrand, thisTopGenSeq)
            if(thisFlag > 0):
                countDifferent += 1
                I_probe = 2
            else:
                countOK += 1
                I_probe = 1
#                 print thisIlmnID + " changed"
        else:
            countOmitted += 1
            if(thisDict2.has_key(thisSnpName)):
                thisFlag2 = compareData(thisDict2, thisSnpName, thisIlmnID, thisIlmnStrand, thisTopGenSeq)
                if(thisFlag2 > 0):
                    countDiffProbe += 1
                    I_probe = 3
            
#             print "SNP is omitted " + thisIlmnID
            # should also check if the same SNP name got a different probe design
        
        outFH.write(thisSnpName + '\t' + str(I_probe) + '\n')
        
print "Different: " + str(countDifferent)
print "Omitted: " + str(countOmitted)
print "No Change: " + str(countOK)
print "Same SnpName, Different Probe: " + str(countDiffProbe)

thisFile2FH.close()
outFH.close()