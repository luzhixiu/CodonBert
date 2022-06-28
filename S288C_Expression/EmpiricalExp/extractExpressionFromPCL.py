import pandas as pd

inputFileName="GSE96997.mapped.pcl"
df=pd.read_csv(inputFileName,sep="\t")



df["GeneName"] = df["ID_REF"].copy()
df["EXP"]  = df['cdc14-3_0min'].copy()
print(df.columns)
df=df[1:]
df[['GeneName','EXP']].to_csv("exp.csv",index=False)