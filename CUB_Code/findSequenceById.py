#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 09:48:45 2019

@author: lu
"""



#idTypeOption is "locus_tag or gene or Gn" *Gn is the keyword for wormbase ID, "raw" means return the whole header
def findSequenceByID(inputFile,idType="locus_tag"):
    print ("Selected id Type: %s"%(idType))
    geneDict=dict()
    from Bio import SeqIO
    records=SeqIO.parse(inputFile, "fasta")
    cnt=0
    mySum=0
    for record in records:
        mySum+=1
        header=str(record.description)
        if idType=="raw":
            geneDict[header]=str(record.seq)
        else:
            startTargetIndex=header.find(str(idType))
            if startTargetIndex<0:
#                print ("couldn't find the target idType")
                cnt+=1
                continue
            startIndex=startTargetIndex+len(idType)+1
            idName=""
            charIndex=startIndex
            while not (header[charIndex]=="]" or header[charIndex]==","):
                idName+=header[charIndex]
                charIndex+=1
            if idName not in geneDict:
                geneDict[idName]=str(record.seq)
    print ("There are %s entries NOT found out of %s"%(cnt,mySum))
    print ("%s distinct record in %s entries"%(len(geneDict),mySum))
    return geneDict



#findSequenceByID("atha.fasta")
    



#for record in SeqIO.parse("atha.fasta", "fasta"):
#    print(record)
