
import os 
import pandas as pd
import subprocess

output = subprocess.check_output('zgrep  -v "^##" Upasana_X2_merged_kggseq-filtered.flt.vcf.gz |grep "#"', shell=True)
output=output.decode("utf-8").replace("\n","").split("\t")




df=pd.read_csv("Upasana_X2_merged_kggseq-filtered.hg38_multiannoDesiredColumns.tsv",sep="\t")
COLS=[x for x in df.columns if "Unn" in x]
Colreplace=dict(zip(COLS, output))
df.rename(columns=Colreplace,inplace=True)


samplenames=list(df.loc[:, 'FORMAT':].columns[1:])
Colreplace=dict(zip(COLS, output))
df.rename(columns=Colreplace,inplace=True)


df2=df[samplenames]
for col in samplenames:
      df2[col]=df2[col].str.split(":",expand=True)[0]


List=[]
for col in samplenames:
      List=List+list(df2[col].unique())
print(set(List))

df2.replace('0|0','0/0',inplace=True)
df2.replace('1|1','1/1',inplace=True)
df2.replace('0|1','0/1',inplace=True)
df2.replace('1|0','0/1',inplace=True)
df2.replace('1/0','0/1',inplace=True)

List=[]
for col in samplenames:
      List=List+list(df2[col].unique())
print(set(List))



Casecount=df2[['X21','X26']].apply(lambda x: x.value_counts(), axis=1)
Casecount.rename(columns={".":"CaseMissi","0/0":"CaseWT","0/1":"CaseHET","1/1":"CaseHOMO"},inplace=True)


Contcount=df2[['X22', 'X23', 'X24', 'X25']].apply(lambda x: x.value_counts(), axis=1)
Contcount.rename(columns={".":"ContMissi","0/0":"ContWT","0/1":"ContHET","1/1":"ContHOMO"},inplace=True)

Case_control_count=pd.concat([Casecount,Contcount],axis=1)
Results=pd.concat([df,Case_control_count],axis=1)

Results.to_csv("Upasana_X2_merged_kggseq-filtered.hg38_multiannoDesiredColumns_GT_Count.csv",index=None)





