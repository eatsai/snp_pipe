#! /usr/bin/python
import sys


class SnpInfo:
    ''' Python class to parse a .bpm.csv file '''
    def __init__(self, alleleFile):
        self.names = []
        self.chr = []
        self.pos = []
        self.refallele = []
        self.altallele = []
        self.strand = []

        alleleFH = open(alleleFile, 'r')
        for thisLine in alleleFH.readlines():
            thisLine = thisLine.rstrip()
            (thisSnpName, thisChr, thisPos, thisRefAllele, thisAltAllele, thisStrand) = thisLine.split("\t")
            self.names.append(thisSnpName)
            self.chr.append(thisChr)
            self.pos.append(thisPos)
            self.refallele.append(thisRefAllele)
            self.altallele.append(thisAltAllele)
            self.strand.append(thisStrand)
            

            

            
            
            
