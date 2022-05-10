#This is a functional library code that convert a fasta file into rank sequences. psudo steps are as below:
#Read in a fasta file downloaded from genbank, use geneId or locus tag to extract the sequences
#Calculate the RSCU for each codn
#Repalce the codons with the RSCU Rank
import sys; sys.path.append('../CUB_Code')
import findSequenceById as FSBID
from CAI import RSCU
import scipy.stats as ss
import CodonLibraries as CL
inputFile="../Data/s288c.fasta"

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

def process(inputF,outputF,tag="", write = False):
    geneDict=FSBID.findSequenceByID(inputF)
    keyList=[]
    seqList=[]
    cnt=0
    for key in geneDict:
        seq=geneDict[key]
        if len(seq)%3==0:
            keyList.append(key)
            seqList.append(seq)
        else:
            cnt+=1
    print(cnt," sequences are not length of three")
    rscu=RSCU(seqList)
    rscu_rank=convertRSCUtoRanks(rscu)
    sentenceList = ['ranks, amino']
    for seq in seqList:
        codonList=CL.loadSequence(seq)
        #remove the first and last five codons:
        codonList = codonList[5:]
        codonList = codonList[:-5]
        # print(codonList)
        # exit()
        try:
            codonRankList=[rscu_rank[codon] for codon in codonList]
            sentence=""

            for rank in codonRankList:
                sentence += str(rank) + " "

            sentence += ', '

            for i in range(len(codonList)):
                add = " " if i < len(codonList) - 1 else ""
                sentence += codonTable[codonList[i]] + add

            sentence += ','


            sentenceList.append(sentence)
        except :
            print("one error on ",seq)

    # print(sentenceList[:2])
    # exit()

    if write:
        output_file = open(outputF, 'w')
        for sentence in sentenceList:
            output_file.write(sentence + '\n')
        output_file.close()




def convertRSCUtoRanks(rscu):
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
        'TYR': ['TAT', 'TAC']}#'CYS': ['TGT', 'TGC']
    rscu_rank=dict()
    for aa in synonymousCodons:
        codonList=synonymousCodons[aa]
        rscuList=[rscu[codon] for codon in codonList]
        rankList=ss.rankdata([-1*x for x in rscuList])
        for codon,rank in zip(codonList,rankList):
            rscu_rank[codon]=int(rank)
    return rscu_rank


def processDir(inputDir, outputDir):
    import os
    files=os.listdir(inputDir)
    for f in files:
        fpath=os.path.join(inputDir,f)
        print('In:', fpath, '\n')
        outputF = os.path.join(outputDir, f + ".CodonRank.txt")
        print('Out:', outputF, '\n')
        process(fpath,outputF, write = True)

processDir("../Data/fnaCollection", outputDir = '../Data/fna_rankC')