#This is a functional library code that convert a fasta file into rank sequences. psudo steps are as below:
#Read in a fasta file downloaded from genbank, use geneId or locus tag to extract the sequences
#Calculate the RSCU for each codn
#Repalce the codons with the RSCU Rank
import findSequenceById as FSBID
from CAI import RSCU
import scipy.stats as ss
import CodonLibraries as CL
inputFile="../Data/s288c.fasta"






def process(inputF,outputF,tag=""):
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
    sentenceList = []
    for seq in seqList:
        codonList=CL.loadSequence(seq)
        #remove the first and last five codons:
        codonList = codonList[5:]
        codonList = codonList[:-5]
        try:
            codonRankList=[rscu_rank[codon] for codon in codonList]
            sentence=""
            for rank in codonRankList:
                sentence+=str(rank)+" "
            sentenceList.append(sentence)
        except :
            print("one error on ",seq)

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


def processDir(inputDir):
    import os
    files=os.listdir(inputDir)
    for f in files:
        fpath=os.path.join(inputDir,f)
        print(fpath)
        outputF=fpath+".CodonRank.txt"
        process(fpath,outputF)

