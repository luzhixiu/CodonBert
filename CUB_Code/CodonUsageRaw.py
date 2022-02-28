#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:36:56 2021

@author: lu
"""
import findSequenceById as FSB
from CAI import RSCU
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import matplotlib.pyplot as plt

targetFastaFile="/home/lu/AcrossTissue/Fastas/c_elegan.fasta"


codonTable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def loadSequence(sequence):
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    actualCodonList=[]
    started=False
    for codon in codonList:
        if codon in stopCodonList:
            break
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList

def RSCU_To_List(rscu):
    resultList=[]
    for key in rscu:
        resultList.append(rscu[key])
    return resultList

def testCorelation(x,y,corelationFunction):
    import scipy.stats as ss
    if "pearson" in corelationFunction:
#        print ("p value %f"%stats.pearsonr(x, y)[1])
        return  ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x,y,nan_policy="omit")[0]


geneDict=FSB.findSequenceByID(targetFastaFile,idType="raw")

print("Adult Set:")
seqList=[]
for gene in geneDict:
    seq=geneDict[gene]
    seqList.append(seq)

adultRSCU= RSCU(seqList)   
correlationList=[]
for gene in geneDict:
    seq=geneDict[gene]
    rscu=RSCU([seq])
    gene_profile=RSCU_To_List(rscu)
    correlationList.append(testCorelation(gene_profile,RSCU_To_List(adultRSCU),"spearman"))
    
#print(correlationList)
plt.figure()
plt.xlabel("Correlation")
plt.xlim(0,1)
plt.hist(correlationList,bins=20)
plt.show()

wholeGenomeFasta="/home/lu/AcrossTissue/Fastas/c_elegan.fasta"
wholeGenomeDict=FSB.findSequenceByID(wholeGenomeFasta,idType="Gn")

wholeGenomeSeqList=[]
for key in wholeGenomeDict:
    seq=wholeGenomeDict[key]
    wholeGenomeSeqList.append(seq)
    
import random

randomGenes=random.sample(wholeGenomeSeqList,1000)

print("Random C elegan Genes")
CorrelationList=[]
for gene in randomGenes:
    try:
        rscu=RSCU([gene])
    except:
        continue
    gene_profile=RSCU_To_List(rscu)
    correlationList.append(testCorelation(gene_profile,RSCU_To_List(adultRSCU),"spearman"))

plt.figure()
plt.xlim(0,1)
plt.xlabel("Correlation")
plt.hist(correlationList,bins=20)
plt.show()
   
def testCorelation(x,y,corelationFunction):
    import scipy.stats as ss
    if "pearson" in corelationFunction:
#        print ("p value %f"%stats.pearsonr(x, y)[1])
        return  ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x,y,nan_policy="omit")[0]
    
testCorelation()
    

    