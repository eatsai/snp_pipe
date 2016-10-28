#! /usr/bin/python
import sys

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

# The Illumina provided Code was provided as-is and with no warranty as to performance and no warranty against it infringing any other party's intellectual property rights.

# Modified 4/6/2015 to work with the 850K BPM.csv by Jeremy Peirce

BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}

class BPM:
    ''' Python class to parse a .bpm.csv file '''
    def __init__(self, bpmFile):
        self.ilmnid = []
        self.names = []
        self.chr = []
        self.pos = []
        self.ilmnStrand = []
        self.normID = []
        self.A = []
        self.B = []
        self.sourceSeq = []
        self.topSeq = []

        with open(bpmFile, 'r') as bpmFilename:
            for _ in xrange(8):  #skips the header range
                next(bpmFilename)
            for line in bpmFilename:
                line = line.replace("\n", "")
                line = line.replace("\r", "")

                fields = line.split(",")

                if line.find('[Controls]') != -1:
#                    print 'found controls'
                    break


                else:
                    self.ilmnid.append(fields[0])   #IlmnID
                    self.names.append(fields[1])  #Name
                    self.chr.append(fields[9])     #Chr
                    self.pos.append(fields[10])    #MapInfo
                    self.ilmnStrand.append(fields[2])
                    alleles = fields[3].replace("[", "") #SNP
                    alleles = alleles.replace("]", "")
                    alleles = alleles.split("/")
                    self.A.append(alleles[0]) # allele A
                    self.B.append(alleles[1]) # allele B
#                     thisIlmnS = fields[2]
                    thisSourceSeq = fields[16]
                    self.sourceSeq.append(thisSourceSeq)
                    thisTopSeq = fields[17]
                    self.topSeq.append(thisTopSeq)
                    
                    
                    ### TODO: delete this conditional statement
#                     if(thisIlmnS == thisSourceS):
#                         self.FwdA.append(alleles[0]) # allele A
#                         self.FwdB.append(alleles[1]) # allele B
#                     elif(thisIlmnS != thisSourceS): # flip strand
#                         self.FwdA.append(BaseComp[alleles[0]]) # allele A
#                         self.FwdB.append(BaseComp[alleles[1]]) # allele A
#                     else:
#                         print "Something went wrong"
#                         sys.exit()
                    self.normID.append(int(fields[18])) # normalization ID for that snp (BeadSetID in newer BPM?)

