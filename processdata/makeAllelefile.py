#! /usr/bin/python

import pysam, re, sys
from BPM import *
from Strand import *
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
    return ('0','0')

def reverseSeq(thisSeq):
    thisRevSeq = ''
    if(thisSeq == '-'):     # TODO: make sure indels don't cause off by one errors later on for VCFs
        return thisSeq
    for thisBase in reversed(thisSeq):
        thisRevSeq += BaseComp[thisBase]
    return (thisRevSeq)

### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-B","--bpm",
                    type=str,
                    required=False,
                    dest="bpm",
                    help="csv-version of bpm file")

parser.add_argument("-S","--strandlist",
                    type=str,
                    required=False,
                    dest="strandlist",
                    help="space-delimited strandlist file (ilmnid, chr, pos, ref allele, strand, I_nonuniq, I_addpos, I_diffpos, I_nomap)")

parser.add_argument("-E","--excludesnplist",
                    type=str,
                    required=False,
                    dest="excludesnplist",
                    help="a list of all the snpnames to exclude (e.g. the SNPs that were zeroed from manual clustering)")

parser.add_argument("-R","--refgenome",
                    type=str,
                    default="/Volumes/pcpgm/et85/share/reference/hg19/Homo_sapiens_assembly19.fasta",
                    dest="refgenome",
                    help="path to reference FASTA file")

# parser.add_argument("-I","--includesnplist",
#                     action='store_true',
#                     required=False,
#                     dest="includesnplist",
#                     help="a list of all the snpnames to include")


args = parser.parse_args()

if(args.bpm):
    bpm = BPM(args.bpm)
if(args.strandlist):
    strandData = Strand(args.strandlist)

if(not (args.bpm and args.strandlist)):
    print "Please make sure to include bpm and strandlist files"

if(args.excludesnplist):
    I_exclude = 1
    excludeSnps = {}
    excludeFH = open(args.excludesnplist,"r")
    for thisLine in excludeFH.readlines():
        thisLine = thisLine.rstrip()
        excludeSnps[thisLine] = 1
else:
    I_exclude = 0

    
### Initialize BPM file class
BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}
numSNPs = len(bpm.names)
# reference_fasta = '/Volumes/pcpgm/share/reference/hg19/Homo_sapiens_assembly19.fasta'
ref_fasta = args.refgenome
pysam_fasta = pysam.FastaFile(ref_fasta)


### Make MAP file
for i in range(numSNPs):
    # assume SNP in bpm file is in the asme order as strand data
    # NB: the new formatting of IlmnID versus SnpName has diverged enough for the check not to work (20160127)
    
    if(bpm.ilmnid[i] == strandData.ilmnid[i]):
        thisSnp = bpm.names[i]
        thisChr = strandData.chr[i]
        thisPos = strandData.pos[i]
    else:
        print "bpm.csv and strandlist files do not appear to be in the same order."
        print bpm.ilmnid[i] + " " + strandData.names[i]
        sys.exit()
    
    if(args.strandlist and strandData.strand[i] == 'u'):
#         continue
        # printing output even for unmapped SNPs
        out = [thisSnp, '0', '0', 'u', bpm.A[i] + ',' + bpm.B[i], 'u'] # zero out chr pos, and put unknown ref/alt; decided to put the alternative TOP SNPs in the ALT allele column
        if(I_exclude):
            if(excludeSnps.has_key(out[0])):
                pass
            else:
                print "\t".join(out)
        else:
            print "\t".join(out)
    elif(args.strandlist and thisPos == '0'):         # there's a new class of data
        out = [thisSnp, '0', '0', 'u', bpm.A[i] + ',' + bpm.B[i], 'u'] # zero out chr pos, and put unknown ref/alt; decided to put the alternative TOP SNPs in the ALT allele column
        if(I_exclude):
            if(excludeSnps.has_key(out[0])):
                pass
            else:
                print "\t".join(out)
        else:
            print "\t".join(out)
        

#         print strandData.strand
    else:   # only print things that have an associated strand        
        (thisA, thisB) = getAlleles(bpm.topSeq[i])
        if(strandData.allele[i] == thisA):
            (tempRef, tempAlt) = getFwdRefAlt(thisA, thisB, strandData.strand[i])
        elif(strandData.allele[i] == thisB):
            (tempRef, tempAlt) = getFwdRefAlt(thisB, thisA, strandData.strand[i])
#         elif(strandData.allele[i] == BaseComp[thisA]): # last two elif should not exist if we BLAT on SourceSeq and match on SourceSeq
#             (thisFwdRef, thisFwdAlt) = getFwdRefAlt(BaseComp[thisA], BaseComp[thisB], strandData.strand[i])
#         elif(strandData.allele[i] == BaseComp[thisB]):
#             (thisFwdRef, thisFwdAlt) = getFwdRefAlt(BaseComp[thisB], BaseComp[thisA], strandData.strand[i])
        else:   # error
            print "Something went wrong..."
        
        
        # now that we got what we think are the ref/alt alleles...
        if(tempRef != '-'):
            if(tempAlt != '-'):
                # dealing with pure SNVs here
                if(thisChr == 'XY'):
                    thisFwdRef = pysam_fasta.fetch('X', int(thisPos) - 1, int(thisPos) - 1 + len(tempRef))
                else:
                    thisFwdRef = pysam_fasta.fetch(thisChr, int(thisPos) - 1, int(thisPos) - 1 + len(tempRef))
                
                if(thisFwdRef == tempRef):
                    thisFwdAlt = tempAlt
                else:   # uh oh, this gets complicated
                    if(thisFwdRef == tempAlt):
                        # not sure why this happens (e.g. "200610-10")
                        # maybe it's because chrMT is unreliable?
                        # flip ref/alt
                        thisFwdRef = tempAlt
                        thisFwdAlt = tempRef
                        
#                         print '\t'.join([thisSnp, thisChr, thisPos, thisFwdRef, thisFwdAlt])
                        
                    else:   #### ERRORS WITH MANIFEST FILE
#                         thisFwdAlt = tempRef + ',' + tempAlt
                        if((tempRef == 'A' and tempAlt == 'T') or (tempRef =='T' and tempAlt == 'A')
                           or (tempRef == 'C' and tempAlt == 'G') or (tempRef == 'G' and tempAlt == 'C')):
                            # Zero out the coordinates. These SNPs are unreliable!
                            thisChr = '0'
                            thisPos = '0'
                            thisFwdRef = 'u'
                            thisFwdAlt = 'u'
                        else:   # TODO: We should flag this, because we only GUESS that the SNP was flipped, but it could be half-flipped or something
                            # we should know how to flip this
                            # 1. which snp is complementary to the forward ref
                            # 2. flip the alternative allele
                            if(thisFwdRef == BaseComp[tempRef]):
                                thisFwdAlt = BaseComp[tempAlt]
                            elif(thisFwdRef == BaseComp[tempAlt]):
                                thisFwdAlt = BaseComp[tempRef]
                            else:
                                sys.exit("Something went wrong. Cannot find proper ref/alt pairing.")
                            

            else:
                # dealing with deletions
                if(thisChr == 'XY'):
                    thisFwdRef = pysam_fasta.fetch('X', int(thisPos) - 2, int(thisPos) - 1 + len(tempRef))
                else:
                    thisFwdRef = pysam_fasta.fetch(thisChr, int(thisPos) - 2, int(thisPos) - 1 + len(tempRef))

                thisFwdAlt = thisFwdRef[0]
                thisPos = str(int(thisPos) - 1)
            
        else:
            # dealing with insertions
            if(thisChr == 'XY'):
                thisFwdRef = pysam_fasta.fetch('X', int(thisPos) - 1, int(thisPos))
            else:
                thisFwdRef = pysam_fasta.fetch(thisChr, int(thisPos) - 1, int(thisPos))
            thisFwdAlt = thisFwdRef + tempAlt
            # position remains unchanged for insertion
        
        
        
        
        ## printing output        
        out = [thisSnp, thisChr, thisPos, thisFwdRef, thisFwdAlt, strandData.strand[i]]
        if(I_exclude):
            if(excludeSnps.has_key(out[0])):
                pass
            else:
                print "\t".join(out)
        else:
            print "\t".join(out)





