#! /usr/bin/python

from BPM import *
from Strand import *
from SnpInfo import *
from argparse import ArgumentParser


### Parse Inputs from Command Line
parser = ArgumentParser()

parser.add_argument("-B","--bpm",
                    type=str,
                    required=True,
                    dest="bpm",
                    default='/Users/et85/testing/gtc/MEGA_Consortium_15063755_B1.csv',
                    help="csv-version of bpm file")

parser.add_argument("-S","--strandlist",
                    type=str,
                    required=False,
                    dest="strandlist",
                    help="space-delimited strandlist file (ilmnid, chr, pos, ref allele, strand, I_nonuniq, I_addpos, I_diffpos, I_nomap)")

parser.add_argument("-A","--allelefile",
                    type=str,
                    required=False,
                    dest="allelefile",
                    help="tab-delimited allele file (ilmnid, chr, pos, refAllele, altAllele(s))")

parser.add_argument("-M","--mapped",
                    action='store_true',
                    required=False,
                    dest="I_mapped",
                    help="A flag to only return values that have a strand")

args = parser.parse_args()

if(args.bpm):
    bpm = BPM(args.bpm)

if(args.strandlist):
    strandData = Strand(args.strandlist)
    if(args.allelefile):
        sys.exit("Please only specify strandlist [-S] or allelefile [-A]")
elif(args.allelefile):
    strandData = SnpInfo(args.allelefile)
else:
    sys.exit("Please only specify at least one: strandlist [-S] or allelefile [-A]")


if(args.I_mapped):
    if(args.strandlist):
        I_writeMapped = 1
    else:
        print "No strandlist provided. Option \"-M\" is not used. Writing MAP file instead."
        I_writeMapped = 0
        
else:   # print a map file if you haven't given a strandlist
    I_writeMapped = 0

### Initialize BPM file class
BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}
numSNPs = len(bpm.names)

### Make MAP file
### NB: If strandlist is provided, must use chr:pos from strandlist (b/c it could have been updated from the bpm.csv positions
for i in range(numSNPs):
    thisSnp = bpm.names[i]
    # check to make sure the order of the files is correct:
    if(bpm.names[i] == strandData.names[i]):
        thisChr = strandData.chr[i]
        thisPos = strandData.pos[i]
    else:
        print "bpm.csv and strandlist files do not appear to be in the same order."
        sys.exit()
    
    if(I_writeMapped):
        if(strandData.strand[i] == 'u'):
            continue
        else:
            out = [thisChr, thisSnp, "0", thisPos]
    else:
        out = [thisChr, thisSnp, "0", thisPos]
    print "\t".join(out)
