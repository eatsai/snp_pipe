#! /usr/bin/python
import sys


BaseComp = {'A':'T','C':'G','G':'C','T':'A','D':'D','I':'I'}

class Strand:
    ''' Python class to parse a .bpm.csv file '''
    def __init__(self, strandFile):
        self.ilmnid = []
#         self.names = []
        self.chr = []
        self.pos = []
        self.allele = []
        self.strand = []
        self.I_nonuniq = []
        self.I_addpos = []
        self.I_diffpos = []
        self.I_nomap = []

        strandFH = open(strandFile, 'r')
        for thisLine in strandFH.readlines():
            thisLine = thisLine.rstrip()
            (thisIlmnid, thisChr, thisPos, thisAllele, thisStrand, thisI_nonuniq, thisI_addpos, thisI_diffpos, thisI_nomap) = thisLine.split(" ")
            self.ilmnid.append(thisIlmnid)
            tempList = thisIlmnid.split("-")
#             self.names.append("-".join(tempList[0:len(tempList)-1]))
            self.chr.append(thisChr)
            self.pos.append(thisPos)
            self.allele.append(thisAllele)
            self.strand.append(thisStrand)
            self.I_nonuniq.append(thisI_nonuniq)
            self.I_addpos.append(thisI_addpos)
            self.I_diffpos.append(thisI_diffpos)
            self.I_nomap.append(thisI_nomap)
        
            

            
            
            
