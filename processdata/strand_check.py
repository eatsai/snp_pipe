#! /usr/bin/python

# Purpose:  Modification of process_forward_strand.py so that it runs well my psl scripts
# Author:   et85, etsai@bwh.harvard.edu

import collections, re, sys
from BPMjmod import *
#from SAM import *  #you can import from the "SAM" or "BLAT" modules
from BLAT import *
from optparse import OptionParser

####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def findRefAllele(thisMatchPosList):
#     minMismatch = ['', 100]
    score = ['', 0]
    for thisMatch in thisMatchPosList:

#         if(minMismatch[1] > thisMatch[3]):
#             minMismatch = [thisMatch[10], thisMatch[3]]
        if(score[1] < thisMatch[4]):
            score = [thisMatch[10], thisMatch[4]]
    
    return score[0]
    
    # only depends on score, because indels have 0 mismatch for both alleles.
#     if(minMismatch[0] != score[0]):
#         print "Something went wrong...",
#         print thisMatchPosList
#     else:
#         return score[0]

def matchSeqToBpm(thisPerfectEntries, blatEntries, thisBpmChr, thisBpmPos, thisLeftOffset, thisRightOffset):
    I_chr_value = 'u'
    I_pos_value = 'u'
    thisAllele = 'u'
    bestBlatStr = 'u'
    I_nonuniq_value = 0
    I_addpos_value = 0
    I_diffpos_value = 0
    I_nomap_value = 0
    
    # Make sure there's a sublist to go through
    if(thisPerfectEntries):
        matchPos = thisPerfectEntries[0][1]
        for testEntry in thisPerfectEntries:
            curPos = testEntry[1]
            if(curPos != matchPos):
                I_nonuniq_value = 1
        if(thisBpmChr == '0' and not I_nonuniq_value):      # no chr/pos to match to
            thisMatchPosList =[]
            I_addpos_value = 1
        elif(thisBpmChr == 'XY'):
            thisMatchPosList = [e for e in thisPerfectEntries if (('X' == e[0] or 'Y' == e[0]) and (abs(thisBpmPos - int(e[1])) <= max([thisLeftOffset,thisRightOffset]) + 1))]
    #       I_matchPos = any((('X' == e[0] or 'Y' == e[0]) and thisBpmPos == e[1]) for e in thisPerfectEntries)
        else:
            thisMatchPosList = [e for e in thisPerfectEntries if (thisBpmChr == e[0] and  (abs(thisBpmPos - int(e[1])) <= max([thisLeftOffset,thisRightOffset]) + 1))]
    #       I_matchPos = any((thisBpmChr == e[0] and thisBpmPos == e[1]) for e in thisPerfectEntries)
    else:
        thisMatchPosList = []
    ## If there's only one position for perfectMatch --> set it.
        # maybe do this later in the logic
        
    if(not thisMatchPosList and blatEntries):     # if this position doesn't match
        # @TODO: need to figure out a way to determine that this BLAT entry is the BEST entry (that would be uniq, not nonuniq)
        matchPos = blatEntries[0][1]
        for testEntry in blatEntries:
            curPos = testEntry[1]
            if(curPos != matchPos):
                I_nonuniq_value = 1
        
        if(thisBpmChr == 'XY'):
            thisMatchPosList = [e for e in blatEntries if (('X' == e[0] or 'Y' == e[0]) and (abs(thisBpmPos - int(e[1])) <= max([thisLeftOffset,thisRightOffset]) + 1))]
#                 I_matchPos = any((('X' == e[0] or 'Y' == e[0]) and thisBpmPos == e[1]) for e in blatEntries)
        else:
            thisMatchPosList = [e for e in blatEntries if (thisBpmChr == e[0] and (abs(thisBpmPos - int(e[1])) <= max([thisLeftOffset,thisRightOffset]) + 1))]
#                 I_matchPos = any((thisBpmChr == e[0] and thisBpmPos == e[1]) for e in blatEntries)
    
    if(thisMatchPosList):                     # entry found, put in original chr:pos
        # now that we found the right position, we need to identify the right strand
        
        # Don't trust Illumina BPM file
#         I_chr_value = thisBpmChr
#         I_pos_value = thisBpmPos

        I_chr_value = thisPerfectEntries[0][0]
        if(bestBlatStr == '+'):
            I_pos_value = int(thisPerfectEntries[0][1]) + thisRightOffset + 1
        elif(bestBlatStr == '-'):
            I_pos_value = int(thisPerfectEntries[0][9]) - thisLeftOffset

        thisAllele = findRefAllele(thisMatchPosList)
        bestBlatStr = thisMatchPosList[0][2]
        # write out whatever was matched
    elif(I_addpos_value):
        if(abs(thisPerfectEntries[0][4]-thisPerfectEntries[0][5]) <= 5
               and abs(thisPerfectEntries[0][5]-thisPerfectEntries[0][6]) <= 5
               and thisPerfectEntries[0][7] <= 1 and thisPerfectEntries[0][8] <= 5):
            I_nonuniq_value = 0         # now that we looked through all thisPerfectEntries and blatEntries and found nothing, we can go back to insert that unique perfect entry...
#             thisAllele = thisPerfectEntries[0][10]
            thisAllele = findRefAllele(thisPerfectEntries)
            bestBlatStr = thisPerfectEntries[0][2]
            I_chr_value = thisPerfectEntries[0][0]
            if(bestBlatStr == '+'):
                I_pos_value = int(thisPerfectEntries[0][1]) + thisRightOffset + 1
            elif(bestBlatStr == '-'):
                I_pos_value = int(thisPerfectEntries[0][9]) - thisLeftOffset
            else:
                print "Something went wrong..."
        else:
            I_nonuniq_value = 0
            I_addpos_value = 0
            I_diffpos_value = 0
            I_nomap_value = 1  
    elif(thisPerfectEntries):                                     # cannot find entry, this is either I_addpos or I_diffpos
        # If you didn't find anything different, let's go back to step one (perfectEntries) and figure out what we can salvage here:
        if(len(thisPerfectEntries) == 1):   # Now that we've confirmed there's no matches with this chr:pos, go back and sublist again
            I_nonuniq_value = 0
        if(not I_nonuniq_value):
            if(abs(thisPerfectEntries[0][4]-thisPerfectEntries[0][5]) <= 5
               and abs(thisPerfectEntries[0][5]-thisPerfectEntries[0][6]) <= 5
               and thisPerfectEntries[0][7] <= 1 and thisPerfectEntries[0][8] <= 5):
                I_diffpos_value = 1
#                 thisAllele = thisPerfectEntries[0][10]
                thisAllele = findRefAllele(thisPerfectEntries)
                bestBlatStr = thisPerfectEntries[0][2]
                I_chr_value = thisPerfectEntries[0][0]      # all entries should have same position if not I_nonuniq_value
                if(bestBlatStr == '+'):
                    I_pos_value = int(thisPerfectEntries[0][1]) + thisRightOffset + 1
                elif(bestBlatStr == '-'):
                    I_pos_value = int(thisPerfectEntries[0][9]) - thisLeftOffset
                else:
                    print "Something went wrong..."
            else:
                I_nonuniq_value = 0
                I_addpos_value = 0
                I_diffpos_value = 0
                I_nomap_value = 1
        else:       # chr pos str = 'u u u' like default
            I_addpos_value = 0
            I_diffpos_value = 0
            I_nomap_value = 0
    else:       # nothing matches!
        I_nonuniq_value = 0
        I_addpos_value = 0
        I_diffpos_value = 0
        I_nomap_value = 1
        
    
    return [I_chr_value, I_pos_value, thisAllele, bestBlatStr, I_nonuniq_value, I_addpos_value, I_diffpos_value, I_nomap_value]



####################################################################################
########################## Parse Inputs from Command Line ##########################
####################################################################################
parser = OptionParser()
parser.add_option("-B","--bpm",type="string",
                  dest="bpm",action="store",
                  help="CSV format of BPM file")
parser.add_option("-A","--alignedResults",type="string",dest="blat",action="store",help="BLAT alignments for SNPs")
parser.add_option("-O","--output",type="string",dest="output",action="store",
                  help="specify output format for further alignment work as text, fasta, or doublefasta, or as a strandlist of matched strands (name, strand) for gtc2ped.  Specifying diagnostics will generate summary statistics")
parser.add_option("-M","--forwardmaster",type="string",dest="fwd_master",action="store",help="strandfile to use as base")

(options, args) = parser.parse_args()

if options.blat == None:
    print "specify blat alignment file with -A"

if options.bpm == None:
    print "specify BPM file path with -B"
    sys.exit()

if options.output == None:
    print "specify output format with -O"
    sys.exit

### Initialize BLAT and BPM file classes
bpm = BPM(options.bpm)
blat = BLAT(options.blat)

### Get Number of SNPs and blat alignments
numSNPs = len(bpm.names)
numBlatAlignments = len(blat.names)
sortedStrandDict = collections.OrderedDict()


# Initialize ordered dictionary:
for i in range(numSNPs):
    sortedStrandDict[bpm.names[i]] = [bpm.names[i],'u','u','u','u',0,0,0,0]


### if there is a master strand list file to update, create the master list.
### if not, create a dummy master list
fwd_masterStrandList = []
if options.fwd_master != None:
    with open (options.fwd_master, 'r') as fwd_masterFilename:
        for line in fwd_masterFilename:
            line = line.rstrip()
            thisList = line.split()
            sortedStrandDict[thisList[0]] = thisList
#             fwd_masterStrandList.append(fields)
# else:           # not sure if this is necessary
#     for i in xrange(numSNPs):
#         fwd_masterStrandList.append([bpm.names[i],'u','u','u',0,0,0,0])



if options.output == "diagnostics":
    print 'Diagnostics Mode.  Read BPM and Alignment File'
    print 'There are ', numSNPs, ' SNPs'
    print 'There are ', numBlatAlignments, ' BLAT alignments'
    
### Blat results Dictionary Setup
snpLookup={}

### Initialize dictionary with an empty list for each value to append later
for name in bpm.names:
    snpLookup.setdefault(name, [])

### Set up dictionary with key=SNP name, values will be blat entries in the format
### (chr, pos, strand, mismatches) for each entry in a list
for i in xrange(numBlatAlignments):
    # @TODO: Catch exceptions where snpLookup does not have this key
#     snpLookup[blat.names[i]].append([blat.chr[i], blat.pos[i], blat.strand[i], blat.mismatches[i], blat.score[i], blat.qsize[i], blat.tsize[i]])
    snpLookup[blat.names[i]].append([blat.chr[i], blat.pos[i], blat.strand[i], blat.mismatches[i], 
                                     blat.score[i], blat.qsize[i], blat.tsize[i], blat.gapcount[i], blat.gapbases[i], blat.end[i], blat.allele[i]])


if options.output == "diagnostics":
    print 'created dictionary'

"""For each SNP name in bpm find the snp name in the dictionary of blat results
Check the chromosome, then position in the blat results.
If they both match, position to within 50bp, then write the strand.
If no match check the next blat result.
If no more blat results write the strand as a u """

strandList = []
numForward = 0
numReverse = 0
numUnknown = 0
I_nonuniq = []  #read is non-uniquely mapping equally well to more than one position
I_addpos = [] #read was 0:0 and we have filled in the chr:pos
I_diffpos = [] #read was uniquely mapped to a chr:pos different from that in the csv.bpm
I_nomap = [] #read does not map at all to any position
I_chr = []
I_pos = []

for i in xrange(numSNPs): #numSNPs
    #starting values reset for each SNP
    maxMismatch = 50 #set max mismatches very high (used in determining if matches are equally good)
    maxScore = 0 #set max score very high (used in determining if matches are equally good)
    numMatch = 0
    bestBlatStr = "x" #if the best strand is not set below, it is unknown
    writeFlag = 0
    I_nonuniq_value = 0
    
    thisBpmChr = bpm.chr[i]
    thisLeftOffset = bpm.sourceStrand[i].index('[')
    thisRightOffset = len(bpm.sourceStrand[i]) - bpm.sourceStrand[i].index(']') - 1
    thisBpmPos = int(bpm.pos[i])
    
    # @TODO: Write handling when you cannot find this SNP in the reads returned
    
    if(snpLookup.has_key(bpm.names[i])):
        blatEntries = snpLookup[bpm.names[i]]  # gets the list of blat results for the ith SNP name
        numEntries = len(blatEntries)
        if not blatEntries : # cannot find SNP, it was not aligned or queried, if it doesn't have existing info, assume it was not aligned.
            if(options.fwd_master):
                writeFlag = 0     
            elif(sortedStrandDict[bpm.names[i]][4] == 'u' and sortedStrandDict[bpm.names[i]][5] != 1):
                I_chr_value = 'u'
                I_pos_value = 'u'
                thisAllele = 'u'
                bestBlatStr = 'u'
                I_nonuniq_value = 0
                I_addpos_value = 0
                I_diffpos_value = 0
                I_nomap_value = 1
                
                # @TODO: fix this for fwd/rev strand chr positions
                sortedStrandDict[bpm.names[i]][1] = I_chr_value
                sortedStrandDict[bpm.names[i]][2] = I_pos_value
                sortedStrandDict[bpm.names[i]][3] = thisAllele
                sortedStrandDict[bpm.names[i]][4] = bestBlatStr # makes a list of blat strands one per SNP
                sortedStrandDict[bpm.names[i]][5] = I_nonuniq_value
                sortedStrandDict[bpm.names[i]][6] = I_addpos_value
                sortedStrandDict[bpm.names[i]][7] = I_diffpos_value
                sortedStrandDict[bpm.names[i]][8] = I_nomap_value
                
        else:
            writeFlag = 1
#                 print bpm.names[i],
#                 print blatEntries
    else:   # TODO: fix this error
        print "Something went wrong here..."
#         continue    # @TODO: might want to print an error somewhere
#         print blatEntries
    
    
    if(writeFlag):
        ## BLAT query only returning 0 or 1 mismatches, but can include gaps
        thisPerfectEntries = []
        
        # look for perfect matches first:
        # if target size = qsize = tsize
        for thisEntry in blatEntries:
            if((thisEntry[4] == thisEntry[5] == thisEntry[6]) and (thisEntry[7] == 0) and (thisEntry[8] == 0)): # perfect match when Score = QSize = TSize
                thisPerfectEntries.append(thisEntry)
    #             if(bpm.names[i] == '200610-10-0_B_R_1867864658'):
    #                 print thisEntry
    #         print "len(PM): " + str(len(thisPerfectEntries))
        
        # Check if this read map is unique
        if(len(thisPerfectEntries) > 1):
            (I_chr_value, I_pos_value, thisAllele, bestBlatStr, I_nonuniq_value, I_addpos_value, I_diffpos_value, I_nomap_value) = matchSeqToBpm(thisPerfectEntries, blatEntries, thisBpmChr, thisBpmPos, thisLeftOffset, thisRightOffset)
            
        elif(len(thisPerfectEntries) == 1):
            (I_chr_value, I_pos_value, thisAllele, bestBlatStr, I_nonuniq_value, I_addpos_value, I_diffpos_value, I_nomap_value) = matchSeqToBpm(thisPerfectEntries, blatEntries, thisBpmChr, thisBpmPos, thisLeftOffset, thisRightOffset)
            
        elif(len(thisPerfectEntries) == 0):      # @TODO: there is likely to be some error, write in handling later
            # Look at all the 1 mismatch entries:
            oneMismatchEntries = [e for e in blatEntries if ((abs(e[5] - e[6]) == 0) and (thisEntry[7] < 2) and (thisEntry[8] < 2))]
            if(oneMismatchEntries):     # if oneMismatchEntries is not empty
                (I_chr_value, I_pos_value, thisAllele, bestBlatStr, I_nonuniq_value, I_addpos_value, I_diffpos_value, I_nomap_value) = matchSeqToBpm(oneMismatchEntries, blatEntries, thisBpmChr, thisBpmPos, thisLeftOffset, thisRightOffset)
    
            else:       # @TODO: fix this condition so that it prioritizes the third bets hits (e.g. by tSize,qSize differences)
                if(blatEntries):    # no perfect entries
                    print blatEntries
                    (I_chr_value, I_pos_value, thisAllele, bestBlatStr, I_nonuniq_value, I_addpos_value, I_diffpos_value, I_nomap_value) = matchSeqToBpm([], blatEntries, thisBpmChr, thisBpmPos, thisLeftOffset, thisRightOffset)

    # If there's no matches, see if the previous match was u u u u 1 0 0 0
    
#     if(writeFlag and not (sortedStrandDict[bpm.names[i]][5] == '1' and I_nomap_value == 1)):
#     if(writeFlag and (sortedStrandDict[bpm.names[i]][1] =='u' or sortedStrandDict[bpm.names[i]][4] == 'u')):      # starting writing differences:
        sortedStrandDict[bpm.names[i]][1] = I_chr_value
        sortedStrandDict[bpm.names[i]][2] = I_pos_value
        sortedStrandDict[bpm.names[i]][3] = thisAllele
        sortedStrandDict[bpm.names[i]][4] = bestBlatStr # makes a list of blat strands one per SNP
        sortedStrandDict[bpm.names[i]][5] = I_nonuniq_value
        sortedStrandDict[bpm.names[i]][6] = I_addpos_value
        sortedStrandDict[bpm.names[i]][7] = I_diffpos_value
        sortedStrandDict[bpm.names[i]][8] = I_nomap_value
        
    
    #gather forward/reverse stats
    if bestBlatStr == "+":
        numForward += 1
    elif bestBlatStr == "-":
        numReverse += 1
    elif bestBlatStr == "u":
        numUnknown += 1
    
    #prepare output based on setting for unmatched only
    if options.output == 'text':
        print bpm.names[i] ,bpm.chr[i], bpm.pos[i], 'blat', blatEntries
    elif options.output == 'fasta':
        print '>'+bpm.names[i]+'\n'+bpm.sourceStrand[i]
    elif options.output == 'doublefasta':
        sourceSeq = re.split('\[|\/|\]',bpm.sourceStrand[i]) #splits to 4 element list
        fivePrimeSeq = sourceSeq[0]

        if sourceSeq[1] == '-':
            alleleASeq = ''
        else:
            alleleASeq = sourceSeq[1]

        if sourceSeq[2] == '-':
            alleleBSeq = ''
        else:
            alleleBSeq = sourceSeq[2]

        threePrimeSeq = sourceSeq[3]

        alleleAStrand = fivePrimeSeq + alleleASeq + threePrimeSeq
        alleleBStrand = fivePrimeSeq + alleleBSeq + threePrimeSeq

        print '>'+bpm.names[i]+'\n'+alleleAStrand
        print '>'+bpm.names[i]+'\n'+alleleBStrand


# print sortedStrandDict['1:74999375-A-C-0_T_F_2304205013']
# print sortedStrandDict['kgp10056353-0_B_R_1905717549']
# print sortedStrandDict['mnp_rs729172_AC-138_T_F_2298921119']
# print sortedStrandDict['rs11540636-120_T_R_1514241389']


if options.output == 'strandlist': #prepare strand list output for use in gtc2ped.py
    for thisKey, thisValue in sortedStrandDict.items():
        print " ".join(str(x) for x in thisValue)

#     for i in range(numSNPs):
#         # check to see if bpm.names exists
#         print i
#         if(bpm.names[i]):
#             #for each i prints a name, strand, and  and info fields
#             print bpm.names[i],I_chr[i],I_pos[i],strandList[i],I_nonuniq[i],I_addpos[i],I_diffpos[i],I_nomap[i]
#         else:
#             print "Something wrong with #" + str(i)

if options.output == 'diagnostics':  
    print 'Forward   : ', numForward
    print 'Reverse   : ', numReverse
    print 'Unknown: ', numUnknown
    print 'Non-Unique : ', I_nonuniq
    print 'Added strand for 0/0 SNP:  ', sum(I_addpos)
    print 'Different Position than BPM:  ', sum(I_diffpos)
    print 'Not aligned at all:  ', sum(I_nomap)
    print I_nonuniq, len(I_addpos), len(I_diffpos), len(I_nomap)
  


