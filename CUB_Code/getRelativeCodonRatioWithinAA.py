import matplotlib.pyplot as plt
# %matplotlib inline

CodonMap = {
'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,'ATA': 0, 'ATG': 0,
'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
'TAA': 0, 'TAG': 0,'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,   'TCG': 0,
'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,'ACT': 0, 'ACC': 0,
'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
'TGT': 0, 'TGC': 0,'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
}
startCodon = "ATG"
stopCodonList = ["TAG", "TAA", "TGA"]

def is_str(v):
    return type(v) is str

def loadSequence(sequence):
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    if not len(sequence)%3==0:
        # print ("Warning! Sequence length not a multiple of 3, will return nothing for this gene" )
        return
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    if not ((codonList[0]) == startCodon and codonList[-1] in stopCodonList):
        # print("Warning! Sequence start or end in invalid codons, will return nothing for this gene")
        return
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

#the a list of sequences string containing all the coding genes should be included, however, single gene analysis can be used in this function as well
#return a map of codon counts
def getCodonCount(seqList):
    codonMap = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
        'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
        'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
        'TAA': 0, 'TAG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
        'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
        'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ACT': 0, 'ACC': 0,
        'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
        'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
        'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
    }
    freqMap=dict()
    #if the input is parsed in as a single string, make it a list of size 1 for iteration.
    if is_str(seqList):
        tempList=[]
        tempList.append((seqList))
        seqList=tempList
    for seq in seqList:
        codonList=loadSequence(seq)
        if codonList==None:
            continue
        for codon in codonList:
            if codon in codonMap:
                codonMap[codon]+=1
            else:
                codonMap[codon]=1
    return codonMap


def getCodonPairCount(seqList,n_codon):
    codonMap=dict()
    startCodon="ATG"
    global stopCodonList

    #if the input is parsed in as a single string, make it a list of size 1 for iteration.
    if is_str(seqList):
        tempList=[]
        tempList.append((seqList))
        seqList=tempList
    for seq in seqList:
        codonList=loadSequence(seq)
        if codonList==None or len(codonList)<100:
            continue
        codon_pair_list=[]
        for i in range(len(codonList)-n_codon+1):
            codonPair=""
            for j in range(n_codon):
                codonPair+=codonList[i+j]
            codonPairHead=codonPair[:3]
            codonPairTail=codonPair[-3:]
            # print(codonPair,codonPairHead,codonPairTail)
            if codonPairHead in stopCodonList:
                # print("caughtHead",codonPairHead)
                continue
            elif codonPairTail == startCodon:
                # print("caughtTail",codonPairTail)
                continue
            else:
                codon_pair_list.append(codonPair)
        for codon in codon_pair_list:
            if codon in codonMap:
                codonMap[codon]+=1
            else:
                codonMap[codon]=1
    return codonMap




def convertCodonCountToRCU(codonCountMap):
    aaRCUMap=dict()
    synonymousCodons = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC']}
    for aa in synonymousCodons:
        # aaRCUMap[aa]=dict()
        tempDict=dict()
        codonCandidates=synonymousCodons[aa]
        for codon in codonCandidates:
            tempDict[codon]=codonCountMap[codon]
        aaRCUMap[aa]=tempDict
    return aaRCUMap

def plotAAmap():

    return





def test():
    sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
    testResult=getCodonCount(sequence)
    convertCodonCountToRCU(testResult)

def test2():
    import loadCodonPairs as LCP
    geneDict=LCP.findSequenceByID("c_elegan.fasta",idType="gene")
    sequenceList=[]
    for geneName in geneDict:
        sequenceList.append(geneDict[geneName])
    # print(sequenceList)
    for i in range(2,9):
        codonPairMap=getCodonPairCount(sequenceList,i)
        # result=convertCodonCountToRCU(testResult)
        print(codonPairMap)
        print(len(codonPairMap))
        distributionList=[]
        for codonPair in codonPairMap:
            codonCount=codonPairMap[codonPair]
            # if codonCount>10:
            distributionList.append(codonCount)
        plt.figure()
        plt.hist(distributionList,log=True)
        plt.gca().set(title=str(i)+' Codon Pair Frequency Histogram', ylabel='Count(log)',xlabel="Codon Pair Count");
        plt.savefig(str(i)+' Codon Pair Frequency Histogram')
        plt.show()


def test3():
    #look at non-appearing codon pairs:
    meaningfulNon_combo=['TTTTAA', 'TTTTAG', 'TTTTGA', 'TTCTAA', 'TTCTAG', 'TTCTGA', 'TTATAA', 'TTATAG', 'TTATGA', 'TTGTAA', 'TTGTAG', 'TTGTGA', 'CTTTAA', 'CTTTAG', 'CTTTGA', 'CTCTAA', 'CTCTAG', 'CTCTGA', 'CTATAA', 'CTATAG', 'CTATGA', 'CTGTAA', 'CTGTAG', 'CTGTGA', 'ATTTAA', 'ATTTAG', 'ATTTGA', 'ATCTAA', 'ATCTAG', 'ATCTGA', 'ATATAA', 'ATATAG', 'ATATGA', 'ATGTAA', 'ATGTAG', 'ATGTGA', 'GTTTAA', 'GTTTAG', 'GTTTGA', 'GTCTAA', 'GTCTAG', 'GTCTGA', 'GTATAA', 'GTATAG', 'GTATGA', 'GTGTAA', 'GTGTAG', 'GTGTGA', 'TATTAA', 'TATTAG', 'TATTGA', 'TACTAA', 'TACTAG', 'TACTGA', 'CATTAA', 'CATTAG', 'CATTGA', 'CACTAA', 'CACTAG', 'CACTGA', 'CAATAA', 'CAATAG', 'CAATGA', 'CAGTAA', 'CAGTAG', 'CAGTGA', 'AATTAA', 'AATTAG', 'AATTGA', 'AACTAA', 'AACTAG', 'AACTGA', 'AAATAA', 'AAATAG', 'AAATGA', 'AAGTAA', 'AAGTAG', 'AAGTGA', 'GATTAA', 'GATTAG', 'GATTGA', 'GACTAA', 'GACTAG', 'GACTGA', 'GAATAA', 'GAATAG', 'GAATGA', 'GAGTAA', 'GAGTAG', 'GAGTGA', 'TCTTAA', 'TCTTAG', 'TCTTGA', 'TCCTAA', 'TCCTAG', 'TCCTGA', 'TCATAA', 'TCATAG', 'TCATGA', 'TCGTAA', 'TCGTAG', 'TCGTGA', 'CCTTAA', 'CCTTAG', 'CCTTGA', 'CCCTAA', 'CCCTAG', 'CCCTGA', 'CCATAA', 'CCATAG', 'CCATGA', 'CCGTAA', 'CCGTAG', 'CCGTGA', 'ACTTAA', 'ACTTAG', 'ACTTGA', 'ACCTAA', 'ACCTAG', 'ACCTGA', 'ACATAA', 'ACATAG', 'ACATGA', 'ACGTAA', 'ACGTAG', 'ACGTGA', 'GCTTAA', 'GCTTAG', 'GCTTGA', 'GCCTAA', 'GCCTAG', 'GCCTGA', 'GCATAA', 'GCATAG', 'GCATGA', 'GCGTAA', 'GCGTAG', 'GCGTGA', 'TGTTAA', 'TGTTAG', 'TGTTGA', 'TGCTAA', 'TGCTAG', 'TGCTGA', 'TGGTAA', 'TGGTAG', 'TGGTGA', 'CGTTAA', 'CGTTAG', 'CGTTGA', 'CGCTAA', 'CGCTAG', 'CGCTGA', 'CGATAA', 'CGATAG', 'CGATGA', 'CGGTAA', 'CGGTAG', 'CGGTGA', 'AGTTAA', 'AGTTAG', 'AGTTGA', 'AGCTAA', 'AGCTAG', 'AGCTGA', 'AGATAA', 'AGATAG', 'AGATGA', 'AGGTAA', 'AGGTAG', 'AGGTGA', 'GGTTAA', 'GGTTAG', 'GGTTGA', 'GGCTAA', 'GGCTAG', 'GGCTGA', 'GGATAA', 'GGATAG', 'GGATGA', 'GGGTAA', 'GGGTAG', 'GGGTGA']
    codonHeadDict=dict()
    for codonPair in meaningfulNon_combo:
        startCodon = "ATG"
        stopCodonList = ["TAG", "TAA", "TGA"]
        codonPairHead = codonPair[:3]
        codonPairTail = codonPair[-3:]
        if codonPairHead in codonHeadDict:
            codonHeadDict[codonPairHead]+=1
        else:
            codonHeadDict[codonPairHead]=1





test2()


