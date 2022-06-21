
import math
import numpy as np
from scipy import optimize
from scipy.optimize import minimize_scalar
from sklearn.preprocessing import minmax_scale
from statistics import mean




deltaEtaFile="Scer_Selection.csv"
deltaMFile="Scer_Mutation.csv"


deltaEtaFile=open(deltaEtaFile)
lines=deltaEtaFile.readlines()
etaDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        etaDict[splitList[1]]=splitList[2]  

deltaMFile=open(deltaMFile)
lines=deltaMFile.readlines()
mDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        mDict[splitList[1]]=splitList[2]
    




codontable = {
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

inverseTable=  {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
                'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
              'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
              'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
              'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
              'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
              'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}
synonymonDict={}
for key in inverseTable:
    valueList=inverseTable[key]
    for value in valueList:
        synonymonDict[value]=valueList


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



#this method removes sequence that cant be handled by mDict, etaDict and "TGG",
def parseSequence(sequence):
    i=0
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    parsedSequence=""
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            if (codon in mDict and codon in etaDict) or (codon in startCodon) or (codon in stopCodonList) and ("TGG" not in codon):
                parsedSequence+=codon
        i+=3
    return parsedSequence
    

def cutSequence(seq):
    sequence=seq
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
    started=FalsewindowSize
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


def roundList(lst):
    roundedList=[]
    decimalPlaces=6
    for i in lst:
        roundedList.append(round(i,decimalPlaces))
    return roundedList

phiDict=dict()
def method4(codonList):
    global mDict
    global etaDict   
    global phiDict
    phiList=[]

    for codon in codonList:
        if codon in phiDict:
             phiList.append(phiDict[codon])
        else:
            maxprob=0.0
            selectedPhi=0.0
            rangeList=[]
            for i in range(1,101):
                rangeList.append(i/100.0)
            for phi in rangeList:
                if codon in mDict:    
                    deltaM=float(mDict[codon])
                else:
                    deltaM=1.0
                if codon in etaDict:
                    deltaEta=float(etaDict[codon])
                else:
                    deltaEta=1.0
                global synonymonDict
                synoList=synonymonDict[codon]
                divisor=np.exp(-1.00*deltaM-(deltaEta*phi))
                dividant=0
                for syno in synoList:
                    if syno in mDict:
                        deltaM=float(mDict[syno])
                    else:
                        deltaM=1.0
                    if syno in etaDict:
                        deltaEta=float(etaDict[syno])
                    else:
                        deltaEta=1.0
                    tmp=np.exp((-1.00*deltaM)-(deltaEta*phi))
    #                print((-1.00*deltaM)-(deltaEta*phi))
    #                print (tmp)
                    dividant+=tmp
                if not dividant==0:
                    prob=divisor/dividant
                else:
                    prob=0
                
                if prob>maxprob:
                    maxprob=prob
                    selectedPhi=phi
            if not selectedPhi==0:
                phiDict[codon]=selectedPhi
                phiList.append(selectedPhi)
            else:
                phiList.append(0.001)
                print ("found 0")
    return phiList
    
def setWindow(inputList,size):
    windowList=[]
    windowSize=size
    cnt=0
    
    while True:
        cnt+=1
        if cnt+windowSize>len(inputList):
            break
        selectedList=inputList[cnt:cnt+windowSize]
        sum=0
        for i in selectedList:
            sum+=float(i)   
        average=sum/len(selectedList)
        windowList.append(average)
    return windowList        

from scipy.stats.mstats import gmean


def cal_mle_windows(sequence,windowSize=10):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)
    windowList=setWindow(phiList,windowSize)
    return windowList
    



def calPhiForGene(sequence):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)
    avg=gmean(phiList)
    return float(avg)



#translate the results in window

def test():
    sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
    codonList=loadSequence(sequence)
    testResults=method4(codonList)
    print(testResults)


#test()
