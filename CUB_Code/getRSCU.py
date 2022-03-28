
Codons = {
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
SynonymousCodons ={
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
'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 'HIS': ['CAT', 'CAC'],
'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
'TRP': ['TGG'],
'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
'GLU': ['GAG', 'GAA'],
'TYR': ['TAT', 'TAC']}

def getRSCU(MySeq):
    MySeq=MySeq.upper ()
    rscuList=[]
    for i in range(0, len(MySeq), 3):
        codon = MySeq [i:i + 3]
        if codon in Codons:
           Codons[codon] += 1
    for aa in SynonymousCodons:
       total=0.0
       rscu = []
       codons = SynonymousCodons[aa]
       for codon in codons:
           total +=Codons[codon]
       denominator=len(SynonymousCodons[aa])
       if total!=0:
          for codon in codons:
              rscu.append((codon,(round((Codons[codon] * len(codons)) /float (total),3)) /denominator))
       rscuList.append(rscu)
#       print rscuList
    return rscuList

