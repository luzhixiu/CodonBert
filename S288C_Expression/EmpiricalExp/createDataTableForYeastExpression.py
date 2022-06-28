import sys

# setting path
sys.path.append('../../CUB_Code/')
import findSequenceById as FSB
import pandas as pd
import numpy as np

def createBinaryList(lst):
    binaryList=[]
    med=np.median(lst)
    for x in lst:
        if x>=med:
            binaryList.append(1)
        else:
            binaryList.append(0)
    return binaryList

geneDict=FSB.findSequenceByID("../../Data/s288c.fasta")
df1=pd.DataFrame(geneDict.items(),columns=["GeneName","Sequence"])
df2=pd.read_csv("exp.csv")
df=pd.merge(df1,df2,on="GeneName")
binList=createBinaryList(df["EXP"])
df["Label"]=binList
df.to_csv("Full_S288C.datatable.csv")


#Create a samll table where middle exp tiers are removed and provide strong exp signals (higher confidence on labels)
topSelectionRatio=0.3
botSelectionRatio=0.3
topSize=int(topSelectionRatio*len(df))
botSize=int(botSelectionRatio*len(df))

lowIndex=botSize
highIndex=len(df)-topSize

#sort the table by exp, and remove the "middle" hunk
df=df.sort_values("EXP")
df=df.drop(df.index[lowIndex:highIndex])

df.to_csv("S288C.datatable_0.3Top_0.3Bot.csv")



