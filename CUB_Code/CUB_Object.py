#This script takes in a fasta file and return a dictionary codons mapped to their CUB object
import findSequenceById as FSBID
from CAI import RSCU
import scipy.stats as ss
import CodonLibraries as CL
inputFile="../Data/s288c.fasta"


class codonObject:
    def __init__(self,codon="NNN"):
        self.codon=codon
        #Define optional information metrics
        self.aa= None
        self.RSCU= None
        self.RSCU_rank=None

#This object contains a dictionary of codons each mapped to their corresponding codon object.
class CUBObject:
    def __init__(self, inputFile):
        codon_map_to_codonObject=dict()
        geneDict = FSBID.findSequenceByID(inputFile)
        keyList = []
        seqList = []
        cnt = 0
        for key in geneDict:
            seq = geneDict[key]
            if len(seq) % 3 == 0:
                keyList.append(key)
                seqList.append(seq)
            else:
                cnt += 1
        print(cnt, " sequences are not length of three")
        rscu = RSCU(seqList)


        rscu_rank,codon_map_to_AA = convertRSCUtoRanks(rscu)
        for codon in rscu_rank:
            co=codonObject(codon)
            co.aa=codon_map_to_AA[codon]
            co.RSCU=rscu[codon]
            co.RSCU_rank=rscu_rank[codon]
            codon_map_to_codonObject[codon]=co
        self.codon_map_to_codonObject=codon_map_to_codonObject



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
            'TYR': ['TAT', 'TAC']}  # 'CYS': ['TGT', 'TGC']
        rscu_rank = dict()
        for aa in synonymousCodons:
            codonList = synonymousCodons[aa]
            rscuList = [rscu[codon] for codon in codonList]
            rankList = ss.rankdata([-1 * x for x in rscuList])
            for codon, rank in zip(codonList, rankList):
                rscu_rank[codon] = int(rank)
        codon_map_to_AA = defaultdict(list)
        for keys, vals in synonymousCodons.items():
            for val in vals:
                codon_map_to_AA[val].append(keys)
        return rscu_rank, codon_map_to_AA




def process(inputF,tag=""):
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
    print(rscu)
    rscu_rank=convertRSCUtoRanks(rscu)
    print(rscu_rank)
    sentenceList = []
    for seq in seqList:
        codonList=CL.loadSequence(seq)
from collections import defaultdict

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
    res = defaultdict(list)
    for keys, vals in synonymousCodons.items():
        for val in vals:
            res[val]=(keys)
    return rscu_rank,res

# Use Case
# inputFile="../Data/s288c.fasta"
# CUB_Object=CUBObject(inputFile)
# process(inputF=inputFile)
# # print(CUB_Object.codon_map_to_codonObject['TGT'].codon)
# print(CUB_Object.codon_map_to_codonObject['TGT'].RSCU)
# print(CUB_Object.codon_map_to_codonObject['TGT'].RSCU_rank)