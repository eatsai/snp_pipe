#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

# Modified from zCall code by Jackie Goldstein by Jeremy Peirce
# Removed rare variant calling code so this now takes a GTC and BPM file (newer format eg CytoSNP 850K) and outputs PED

import re, sys
from GTC import *
from BPM import *
from Strand import *
from SnpInfo import *
from argparse import ArgumentParser
####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def getAlleles(thisSeq):
    seqList = re.split('\[|\/|\]', thisSeq)
    thisA = seqList[1]
    thisB = seqList[2]
    return(thisA, thisB)
    
def getFwdRefAlt(thisRef, thisAlt, thisStrand):
    if(thisStrand == '+'):
        return (thisRef, thisAlt)
    elif(thisStrand == '-'):
        return (reverseSeq(thisRef), reverseSeq(thisAlt))
#         # 2015-12-23: Don't reverse if ambiguous SNPs, because our logic already took care of this when identifying reference strand
#         if(thisRef + thisAlt == 'CG' or thisRef + thisAlt == 'GC'
#            or thisRef + thisAlt == 'AT' or thisRef + thisAlt == 'TA'):
#             return (thisRef, thisAlt)
#         else:
#             return (reverseSeq(thisRef), reverseSeq(thisAlt))
    return ('0','0')

def getTopAllele(thisInd, thisBpm):
    
    ## what to do if it is indel?
    (thisA, thisB) = ('', '')
    
    if(thisBpm.ilmnStrand[thisInd] == 'PLUS' or thisBpm.ilmnStrand[thisInd] == 'MINUS'):
        thisA = thisBpm.A[thisInd]
        thisB = thisBpm.B[thisInd]
#     elif(thisBpm.ilmnStrand[thisInd] == 'MINUS'):
#         thisA = 'D'
#         thisB = 'I'
#     if(thisBpm.A[thisInd] == 'I' or thisBpm.A[thisInd] == 'D'): # indels aren't affected by TOP/BOT
#         thisA = thisBpm.A[thisInd]
#         thisB = thisBpm.B[thisInd]
    elif(thisBpm.ilmnStrand[thisInd] == 'TOP'):
        thisA = thisBpm.A[thisInd]
        thisB = thisBpm.B[thisInd]
    elif(thisBpm.ilmnStrand[thisInd] == 'BOT'):
        thisA = BaseComp[thisBpm.A[thisInd]]
        thisB = BaseComp[thisBpm.B[thisInd]]
    else:
        print 'Something went wrong, check strand '+thisSnp+':'+thisStrand
        thisA = '0'
        thisB = '0'
    
    return(thisA, thisB)
    

def makeGenderTable(thisGenderTablePath, thisGenderDict):
    thisFH = open(thisGenderTablePath, 'r')
    
    for thisLine in thisFH.readlines():
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split('\t')
        thisGenderDict[thisRead[0]] = str(thisRead[2])      # key: SentrixID, value: BiobankID

def reverseSeq(thisSeq):
    thisRevSeq = ''
    if(thisSeq == '-'):     # TODO: make sure indels don't cause off by one errors later on for VCFs
        return thisSeq
    for thisBase in reversed(thisSeq):
        thisRevSeq += BaseComp[thisBase]
    return (thisRevSeq)

def makeConcordanceTable(idtablepath, thisDict):
    thisFH = open(idtablepath, 'r')
    
    for thisLine in thisFH.readlines()[1:]:      # skip first line (header)
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split('\t')
        thisDict[thisRead[1]] = thisRead[0]      # key: SentrixID, value: BiobankID

####################################################################################
################################   Argument Parser   ###############################
####################################################################################
parser = ArgumentParser()

parser.add_argument("-B","--bpm",
                    type=str,
                    required=True,
                    dest="bpm",
                    default='/Users/et85/testing/gtc/MEGA_Consortium_15063755_B1.csv',
                    help="csv-version of bpm file")

parser.add_argument("-S","--strandlist",
                    type=str,
                    dest="strandlist",
                    help="strandlist file generated from strand_check.py")

parser.add_argument("-A","--allelefile",
                    type=str,
                    dest="allelefile",
                    help="allele file generated from makeAllelefile.py")

parser.add_argument("-G","--gtcfile",
                    type=str,
                    dest="gtc",
                    help="input GTC file to decode")

parser.add_argument("-T","--top",
                    action='store_true',
                    required=False,
                    dest="I_top",
                    help="specify -T or --top if you want to report top allele")

parser.add_argument("-C","--idtable",
                    type=str,
                    dest="idtable",
                    help="specify -C or --idtable if you want to put in the idtable")

### put in gender information as well???
# /data/pcpgm/et85/projects/mega/ped_fwdstr_merge/updatesex12_x4931.list
parser.add_argument("-E","--gendertable",
                    type=str,
                    dest="gendertable",
                    help="specify -E or --gendertable if you want to put in the gendertable")



args = parser.parse_args()


if(not args.bpm):
    print "specify BPM file path with -B"
    sys.exit()


I_allelefile = bool(args.allelefile)
I_strandlist = bool(args.strandlist)
I_top = bool(args.I_top)
I_idtable = bool(args.idtable)
I_gendertable = bool(args.gendertable)

if( (I_allelefile and not (I_strandlist or I_top)) or
    (I_strandlist and not (I_allelefile or I_top)) or
    (I_top and not (I_strandlist or I_allelefile)) ):
    # have an allele to work from
    pass
else:
    sys.exit("Please specify EITHER (1) strandlist file path with -S, (2) allele file path with -A, or (3) write as top allele -T")

if(I_strandlist):
    strandData = Strand(args.strandlist)
elif(I_allelefile):
    alleleData = SnpInfo(args.allelefile)

if(I_idtable):
    idTable = args.idtable
    
if(I_gendertable):
    genderTablePath = args.gendertable


# if(not (bool(args.strandlist) ^ args.I_top)):
#     print "Please specify EITHER (1) strandlist file path with -S, (2) allele file path with -A, or (3) write as top allele -T"
#     sys.exit()

if(not args.gtc):
    print "specify GTC file path with -G"
    sys.exit()

# if(not args.strandlist):
#     I_fwd = 0       # print TOP allele, not FORWARD allele
# else:
#     I_fwd = 1
# #     strandFH = open(args.strandfile, "r")
#     strandData = Strand(args.strandlist)

# print "fwd allele? " + str(options.I_fwd) 

### Parse sample name from input gtc file name
sampleName = args.gtc.split("/")
sampleName = sampleName[len(sampleName) - 1]
sampleName = sampleName.split(".")[0]

### Initialize GTC and BPM file classes
bpm = BPM(args.bpm)    
gtc = GTC(args.gtc, 'gtonly')      # TODO: need to fix bpm.normID

### Get Number of SNPs
if(gtc.numSNPs == len(bpm.names)):
    if(I_strandlist):
        if(gtc.numSNPs == len(strandData.names)):
            numSNPs = gtc.numSNPs
        else:
            print "Inconsistent nSNPs between strandlist and gtc+bpm.csv files"
            sys.exit()
    elif(I_allelefile):
        if(gtc.numSNPs == len(alleleData.names)):
            numSNPs = gtc.numSNPs
        else:
            print "Inconsistent nSNPs between allelefile and gtc+bpm.csv files"
            sys.exit()
    elif(I_top):
        numSNPs = gtc.numSNPs
else:
    print "Inconsistent nSNPS between gtc and bpm.csv files"
    sys.exit()

### Make PED file
BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}

### See if there's a concordance table to use to fix the sample information

# put this in concordance table
IdDict = {}
GenderDict = {}
out = [sampleName, sampleName, "0", "0", "-9", "-9"] #output holder in python list; have no sample information so use "0" and "-9" for mid, pid, gender, case/control status


if(I_idtable):
    makeConcordanceTable(idTable, IdDict)
    if(IdDict.has_key(sampleName)):
        biobankID = IdDict[sampleName]
        out[1] = biobankID          # put in new ID if you have it
    else:
        sys.stderr.write('No value found in correspondence table, SenxtrixID = IID = FID = ' + sampleName + '\n')
    
if(I_gendertable):
    makeGenderTable(genderTablePath, GenderDict)
    out[4] = GenderDict[sampleName]


for i in range(numSNPs):
    if(I_strandlist):
        thisSnp = strandData.names[i]
        thisMappedAllele = strandData.allele[i]     # For a given SNP, thisMappedAllele maps to hg19 better than the other allele
        (thisA, thisB) = getAlleles(bpm.sourceSeq[i])
        thisStrand = strandData.strand[i]
        
        if(args.strandlist and thisStrand == 'u'):    # skip the GT for regions we cannot resolve fwd strand
            continue
        elif(bpm.A[i] == 'D' or bpm.A[i] == 'I'):   # handle indels differently
            A = bpm.A[i]
            B = bpm.B[i]
        elif(thisMappedAllele == thisA):
            (thisFwdRef, thisFwdAlt) = getFwdRefAlt(thisA, thisB, thisStrand)
            A = thisFwdRef
            B = thisFwdAlt
        elif(thisMappedAllele == thisB):
            (thisFwdRef, thisFwdAlt) = getFwdRefAlt(thisB, thisA, thisStrand)
            A = thisFwdAlt
            B = thisFwdRef
        else:   # something went wrong
            print 'Something wrong here, check mapped allele correspondence to strandlist: '+thisSnp+':'+thisStrand+':'+thisMappedAllele
            A = '0'
            B = '0'
    elif(I_allelefile):
        thisSnp = alleleData.names[i]
        ## currently, all [I/D] SNPs will have at least one allele as "-" so we can assume the following.        
        thisFwdRef = alleleData.refallele[i]
        thisFwdAlt = alleleData.altallele[i]
        
        (thisTopA, thisTopB) = getTopAllele(i, bpm)
        
        thisStrand = alleleData.strand[i]
        
        # check to see if there are any weird issues:
        if(thisFwdRef == 'u'):      # This SNP was not mappable or ambigious. just use whatever is designated as TOP allele
            A = thisTopA
            B = thisTopB
        elif(thisFwdRef == '-' or thisFwdAlt == '-'):
            # PLINK encodes everything as a single letter
            # I/D SNPs aren't affected by mapping
            A = thisTopA
            B = thisTopB
        else:   # SNVs
            # 2015-12-23: check for ambiguous SNPs
            if(thisStrand == '+'):
                A = thisTopA
                B = thisTopB
            elif(thisStrand == '-'):
                A = BaseComp[thisTopA]
                B = BaseComp[thisTopB]
            else:
                print 'Something went wrong, check strand '+thisSnp+':'+thisStrand
                A = '0'
                B = '0'
#             if((thisFwdRef == thisBpmA) or (thisFwdRef == thisBpmB)):  # TOP == FWD
#                 A = thisBpmA
#                 B = thisBpmB
#             else:   # TOP == REV
#                 A = BaseComp[thisBpmA]
#                 B = BaseComp[thisBpmB]
        
    elif(I_top):   # mapping TOP allele
        (A, B) = getTopAllele(i, bpm)
#         A = bpm.A[i] # Get what A allele is from bpm.csv file [A,G,T,C]
#         B = bpm.B[i] # Get what B allele is from bpm.csv file [A,G,T,C]
    
    origCall = gtc.genotypes[i] # 0 - "No Call", 1 - AA, 2 - AB, 3 - BB

    if origCall == 0:   ## NC
        out.append("0")
        out.append("0")
    elif origCall == 1: ## AA
        out.append(A)
        out.append(A)
    elif origCall == 2: ## AB
        out.append(A)
        out.append(B)
    elif origCall == 3: ## BB
        out.append(B)
        out.append(B)
    else:
        print "Something wrong here"
        sys.exit()

## Output to std out the new calls in PED format
print " ".join(out)




