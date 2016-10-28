#!/usr/bin/env python
'''
Purpose:       Have a list of keys that are the only things you want to extract from a larger file

Usage:         python ped2vcf.py -P infile.ped -O outfile.vcf -M file.map -S strandinfo.txt (-H vcfheader.txt)
 
Author:        et85, etsai@bwh.harvard.edu

Svn Version:   $Id: ped2vcf.py 2947 2015-12-21 22:34:05Z et85 $

'''

from argparse import ArgumentParser
import sys, time

BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}
####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def loadSnpInfo(thisFile):
    thisDict = {}
    thisFH = open(thisFile, "r")
        
    for thisLine in thisFH.readlines():
        thisLine = thisLine.rstrip()
        (thisSnpName, thisChr, thisPos, thisRef, thisAlt, thisStrand, thisCallRateBin, thisRsid, thisVariantEffect) = thisLine.split('\t')
        thisDict[thisSnpName] = [thisChr, thisPos, thisRef, thisAlt, thisStrand, thisCallRateBin, thisRsid, thisVariantEffect]
    thisFH.close()
    
    return thisDict

def loadMap(thisFile):
    thisList = []
    thisFH = open(thisFile, "r")
    thisInd = 0
    for thisLine in thisFH.readlines():
        thisLine = thisLine.rstrip()
        thisSnpName = thisLine.split('\t')[1]
        thisChr = thisLine.split('\t')[0]
        thisPos = thisLine.split('\t')[3]
        thisList.append([thisSnpName, thisChr, thisPos, thisInd])
        thisInd = thisInd + 1
    thisFH.close()
    
    return thisList

def loadPed(thisLine):
    thisLine = thisLine.rstrip()
    thisRead = thisLine.split(' ')
    thisSampName = thisRead[0] + '-' + thisRead[1]
    thisGender = thisRead[4]
    del thisRead[0:6]   # we get gender from CONSTRAK, doesn't make sense to deliver it back to RPDR    
    return (thisSampName, thisRead)

def calculateAC(thisPedSnp1,thisPedSnp2,thisRefSnp,thisAltSnps):
    thisAltSnpList = thisAltSnps.split(',')
    
    thisAC1 = 0
    
    if(thisPedSnp1 == "0" and thisPedSnp2 == "0"):
        thisAC1 = -1
    
    if(thisAltSnpList[0] == thisPedSnp1):
        thisAC1 = thisAC1 + 1
    if(thisAltSnpList[0] == thisPedSnp2):
        thisAC1 = thisAC1 + 1
    
    if(len(thisAltSnpList) == 1):
        return str(thisAC1)
    
    elif(len(thisAltSnpList) == 2):
        thisAC2 = 0
        
        if(thisPedSnp1 == "0" and thisPedSnp2 == "0"):
            thisAC2 = -1
        
        if(thisAltSnpList[1] == thisPedSnp1):
            thisAC2 = thisAC2 + 1
        if(thisAltSnpList[1] == thisPedSnp2):
            thisAC2 = thisAC2 + 1
        
        thisAC = str(thisAC1) + ',' + str(thisAC2)
        return thisAC

def convertAC_to_GT(thisAC):
    thisACList = thisAC.split(',')
    if(len(thisACList) == 1):
        thisAC1 = int(thisACList[0])
        if(thisAC1== -1):
            return './.'
        elif(thisAC1 == 0):
            return '0/0'
        elif(thisAC1 == 1):
            return '0/1'
        elif(thisAC1 == 2):
            return '1/1'
        
    elif(len(thisACList) == 2):
        
        if(thisACList[0] == '-1' and thisACList[1] == '-1'):
            return './.'
        thisAC1 = int(thisACList[0])
        thisAC2 = int(thisACList[1])
        if(thisAC1 + thisAC2 == 0):   # no alternative alleles
            return '0/0'
        elif(thisAC1 + thisAC2 == 1):
            if(thisAC1 == 1):
                return '0/1'
            elif(thisAC2 == 1):
                return '0/2'
        elif(thisAC1 + thisAC2 == 2): # no ref allele
            if(thisAC1 == 1 and thisAC2 == 1):
                return '1/2'
            elif(thisAC1 == 2):
                return '1/1'
            elif(thisAC2 == 2):
                return '2/2'

        

def printDataset(thisSnpList, thisSnpDict, thisOutFH, thisSnpGenoTable, thisSampList):
    
    thisOutFH.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + thisSampList) + '\n')
    
    for thisInd in range(0,len(thisSnpList)):
        [thisSnpName, thisChr, thisPos, nil1, nil2, thisSnpInd] = thisSnpList[thisInd]
        [thisChr, thisPos, fwdAlleleRef, fwdAlleleAlt, thisStrand, thisCallRateBin, thisRsid, thisVariantEffect] = thisSnpDict[thisSnpName]
        
        thisFwdAlleleAltList = fwdAlleleAlt.split(',')
        (fwdAlleleAlt1, fwdAlleleAlt2) = (None, None)
        
        if(len(thisFwdAlleleAltList) == 1):
            fwdAlleleAlt1 = thisFwdAlleleAltList[0]
        elif(len(thisFwdAlleleAltList) == 2):
            (fwdAlleleAlt1, fwdAlleleAlt2) = thisFwdAlleleAltList
            # Exception for now:
            if(fwdAlleleAlt1 == fwdAlleleAlt2):
                fwdAlleleAlt2 = None   
                fwdAlleleAlt =  fwdAlleleAlt1
                thisVariantEffect = thisVariantEffect.split(',')[0]     # assume the other half is a duplicate
        else:
            print "Something went wrong"
        
        thisInfo = ''    
        
        if(thisRsid != '.'):
            thisInfo = thisInfo + 'RSID=' + thisRsid + ';'
        
        if(thisVariantEffect != '.'):
            thisVariantEffect = thisVariantEffect.replace(' ','_')
            thisInfo = thisInfo + 'VariantEffect=' + thisVariantEffect + ';'
        
        if(thisInfo != ''): # strip the last character (extra semicolon)
            thisInfo = thisInfo[:-1]
        else:
            thisInfo = '.'
        
        
        # writing output
        if(not (thisChr == '0' or thisPos == '0')):
            thisOutFH.write('\t'.join([thisChr, thisPos, thisSnpName, fwdAlleleRef, fwdAlleleAlt, '.', '.', thisInfo, 'GT'] + thisSnpGenoTable[thisSnpInd]) + '\n')
#             print '\t'.join([thisChr, thisPos, thisSnpName, fwdAlleleRef, fwdAlleleAlt, '.', '.', thisInfo, 'GT'] + thisSnpGenoTable[thisSnpInd])
#             print '\t'.join([thisChr, thisPos, thisSnpName, fwdAlleleRef, fwdAlleleAlt]),
#             print thisSnpGenoTable[thisSnpInd]


def processSubject(thisRead, thisSnpList, thisSnpDict, thisOutFH, thisNSamp, thisSnpGenoTable):
    
    for thisInd in range(0,len(thisRead)/2):
        (thisSnpName, thisChr, thisPos, nil1, nil2, thisSnpInd) = thisSnpList[thisInd]
        pedSnp1 = thisRead[thisSnpInd*2]
        pedSnp2 = thisRead[thisSnpInd*2 + 1]
        
        [thisChr, thisPos, fwdAlleleRef, fwdAlleleAlt, thisStrand, thisCallRateBin, thisRsid, thisVariantEffect] = thisSnpDict[thisSnpName]
         
        thisFwdAlleleAltList = fwdAlleleAlt.split(',')
        (fwdAlleleAlt1, fwdAlleleAlt2) = (None, None)
         
        if(len(thisFwdAlleleAltList) == 1):
            fwdAlleleAlt1 = thisFwdAlleleAltList[0]
        elif(len(thisFwdAlleleAltList) == 2):
            (fwdAlleleAlt1, fwdAlleleAlt2) = thisFwdAlleleAltList
            # Exception for now:
            if(fwdAlleleAlt1 == fwdAlleleAlt2):
                fwdAlleleAlt2 = None
                fwdAlleleAlt =  fwdAlleleAlt1
#                 thisVariantEffect = thisVariantEffect.split(',')[0]     # assume the other half is a duplicate
        else:
            print "Something went wrong"
        
        
        if(thisChr == 0):   # skip SNPs that don't have a reported position
            continue
        elif(pedSnp1 == 'I' or pedSnp1 == 'D'):
            
            if(len(fwdAlleleRef) < len(fwdAlleleAlt)):  # Ref = D; Alt = I
                thisAC = calculateAC(pedSnp1,pedSnp2,'D','I')
                thisGT = convertAC_to_GT(thisAC)
                
            elif(len(fwdAlleleRef) > len(fwdAlleleAlt)):    # Ref = I; Alt = D
                thisAC = calculateAC(pedSnp1,pedSnp2,'I','D')
                thisGT = convertAC_to_GT(thisAC)
                
            else:
                print "Not sure how to resolve this."
                print thisSnpName + '\t' + fwdAlleleRef + '\t' + fwdAlleleAlt
                sys.exit()
        
        else:   # Note: SNPs are already in fwd strand orientation
            thisAC = calculateAC(pedSnp1,pedSnp2,fwdAlleleRef,fwdAlleleAlt)
            thisGT = convertAC_to_GT(thisAC)
        
        if(thisNSamp == 1):
            thisSnpGenoTable[thisSnpInd] = [thisGT]
        else:
            thisSnpGenoTable[thisSnpInd].append(thisGT)
        
        
        ## sorted lists uses less memory than hash
#         if(thisNSamp == 1): # must initialize for first sample
#             thisSnpGenoTable = lists(len)
#             thisSnpGenoTable[thisSnpInd] = [thisGT]
#         else:   # append to list of list (e.g. table)
#             thisSnpGenoTable[thisSnpInd].append(thisGT)
        
    return thisSnpGenoTable

####################################################################################
################################   Argument Parser   ###############################
####################################################################################
parser = ArgumentParser()

parser.add_argument("-P","--pedfile",
                    type=str,
                    required=True,
                    dest="pedfile",
                    help="input PED file")

parser.add_argument("-M","--mapfile",
                    type=str,
                    required=True,
                    dest="mapfile",
                    help="input MAP file")

parser.add_argument("-H","--header",
                    type=str,
                    dest="header",
                    help="Header of the VCF file EXCLUDING the line starting with #CHROM")

parser.add_argument("-O","--output",
                    type=str,
                    required=True,
                    dest="outfile",
                    help="output VCF file")

parser.add_argument("-S","--strandinfo",
                    type=str,
                    required=True,
                    dest="strandinfo",
                    help="strand info file is a tab-delimited form of the following: SnpName, Chr, Pos, Ref, Alt, RSID, Annotation")

args = parser.parse_args()

if(args.pedfile):
    pedfile = args.pedfile

if(args.mapfile):
    mapfile = args.mapfile

if(args.header):
    header = args.header
else:
    thisFormDate = time.strftime('%Y%m%d')
    header = '##fileformat=VCFv4.2\n##filedate=' + thisFormDate + '\n##source=PED\n'
    header = header + '##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n'
    header = header + '##INFO=<ID=RSID,Number=.,Type=String,Description=\"dbSNP ID in build 147\">\n'
    header = header + '##INFO=<ID=VariantEffect,Number=.,Type=String,Description=\"Variant annotation (if present) in the following format: '
    header = header + 'GeneSymbol|RefSeqTranscriptID:c.Nomenclature|RefSeqProteinID:p.Nomenclature|VariantEffect.\">\n'
    header = header + '##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=3,length=198022430>\n'
    header = header + '##contig=<ID=4,length=191154276>\n##contig=<ID=5,length=180915260>\n##contig=<ID=6,length=171115067>\n'
    header = header + '##contig=<ID=7,length=159138663>\n##contig=<ID=8,length=146364022>\n##contig=<ID=9,length=141213431>\n'
    header = header + '##contig=<ID=10,length=135534747>\n##contig=<ID=11,length=135006516>\n##contig=<ID=12,length=133851895>\n'
    header = header + '##contig=<ID=13,length=115169878>\n##contig=<ID=14,length=107349540>\n##contig=<ID=15,length=102531392>\n'
    header = header + '##contig=<ID=16,length=90354753>\n##contig=<ID=17,length=81195210>\n##contig=<ID=18,length=78077248>\n'
    header = header + '##contig=<ID=19,length=59128983>\n##contig=<ID=20,length=63025520>\n##contig=<ID=21,length=48129895>\n'
    header = header + '##contig=<ID=22,length=51304566>\n##contig=<ID=X,length=155270560>\n##contig=<ID=Y,length=59373566>\n'
    header = header + '##contig=<ID=MT,length=16569>\n'
    
if(args.outfile):
    outfile = args.outfile

if(args.strandinfo):
    strandinfo = args.strandinfo


####################################################################################
################################    Main Function    ###############################
####################################################################################

nSamp = 0

outFH = open(outfile, "w")
outFH.write(header)
# outFH.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',thisSampName]) + '\n')

# build a hash for strandinfo
SnpDict = loadSnpInfo(strandinfo)

# load map file
snpList = loadMap(mapfile)

sampList = []
snpGenoTable = []

# make new sample list based on SNP dictionary
tempList = []
for thisSnp in snpList:
    [thisSnpName, nil1, nil2, thisSnpInd] = thisSnp
    [thisChr, thisPos, fwdAlleleRef, fwdAlleleAlt, nil1, nil2, nil3, nil4] = SnpDict[thisSnpName]
    tempList.append([thisSnpName, thisChr, thisPos, fwdAlleleRef, fwdAlleleAlt, thisSnpInd])
    
sortedSnpList = sorted(tempList, key = lambda x : (x[1],len(x[2]),x[2],x[3],x[4]))

# load ped file
pedFH = open(pedfile,"r")
for thisLine in pedFH.readlines():
    nSamp = nSamp + 1
    
    (thisSampName, thisRead) = loadPed(thisLine)
    sampList.append(thisSampName)
    if(nSamp == 1):
        snpGenoTable = [None] * (len(thisRead)/2)
    
    snpGenoTable = processSubject(thisRead, sortedSnpList, SnpDict, outFH, nSamp, snpGenoTable)


printDataset(sortedSnpList, SnpDict, outFH, snpGenoTable, sampList)


outFH.close()
pedFH.close()
