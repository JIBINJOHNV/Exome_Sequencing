#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
from difflib import SequenceMatcher
from fuzzywuzzy import fuzz
import numpy as np
import argparse


# In[ ]:



parser=argparse.ArgumentParser(description="It is for generating Manhatton and qq plot from genesis software")
parser.add_argument('-AnnovarFile','--AnnovarFile', help="Files generated from the annovar", required=True)
parser.add_argument('-Phenotype','--Phenotype', help="File with phenotype; should cotain column named Phenotyps", required=False)
args=parser.parse_args()

AnnovarFile=args.AnnovarFile
Phenotype=args.Phenotype


# In[ ]:


df=pd.read_csv(AnnovarFile,sep="\t")
phentype=pd.read_csv(Phenotype)


# In[28]:


#df=pd.read_csv("KT41093424W.hg38_multiannoDesiredColumns.tsv",sep="\t")
#phentype=pd.read_csv("phenotype.txt")


# In[29]:


dfHpo=df[['Gene.refGene','HPO_name.refGene']]
dfHpo['HPO_name.refGene']=dfHpo['HPO_name.refGene'].fillna("No HPO")


# In[30]:


def similar(a, b):
    tem=fuzz.partial_ratio(a.lower(), b.lower()) 
    return tem


# In[22]:



for j in range(df.shape[0]):
    Dict={}
    Score=0
    for i in range(phentype.shape[0]):
        test=similar(phentype.iloc[i,0], 
                     ", ".join(dfHpo.loc[j,"HPO_name.refGene"].split(":")[0].split(";")))
        
        Dict[phentype.iloc[i,0]]=test
        Score=Score+int(test)
    Score=Score/(phentype.shape[0]*100)*100
    df.loc[j,"PhenotpeScoreDict"]=str(Dict)
    df.loc[j,"PhenotpeScor"]=Score

         


# In[23]:


df.to_csv(AnnovarFile[:-4]+"_Phenotype.tsv",sep="\t")


# In[ ]:




