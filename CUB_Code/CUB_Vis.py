import findSequenceById as FSBID
import scipy.stats as ss
import CodonLibraries as CL
import CUB_Object
import matplotlib.pyplot as plt
from matplotlib import markers
# get all possible shapes
all_shapes = markers.MarkerStyle.markers.keys()

inputFile="../Data/s288c.fasta"







CUB_Object=CUB_Object.CUBObject(inputFile)
codon_map_to_codonObject=CUB_Object.codon_map_to_codonObject

for codon in codon_map_to_codonObject:
    codon_object=codon_map_to_codonObject[codon]
    rscu=codon_object.RSCU
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
    sentenceList = []
    import pandas as pd
    import seaborn as sns

    for seq in seqList[:3]:
        codonList = CL.loadSequence(seq)
        aaList = [codon_map_to_codonObject[codon].aa for codon in codonList]
        rscuList = [codon_map_to_codonObject[codon].RSCU for codon in codonList]
        df = pd.DataFrame()
        df["AA"] = aaList
        df["RSCU"] = rscuList
        df["Codon"] = codonList
        df["Index"] = list(range(0, len(codonList)))
        g = sns.lineplot(x='Index', y="RSCU", data=df, hue='AA', marker='o')

        plt.show()

    aa=codon_object.aa
    print(aa,codon)

