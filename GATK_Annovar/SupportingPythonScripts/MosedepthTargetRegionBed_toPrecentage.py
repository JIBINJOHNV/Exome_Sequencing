

import pandas as pd
import glob
import argparse
pd.options.display.float_format = "{:,.2f}".format


# 

parser=argparse.ArgumentParser(description="It is for summarising Mosedepth TargetRegion.thresholds outpt")
parser.add_argument('-Mfolder','--MosedepthOutPutFolder', help="Mosedepth Output Folder", required=True)
args=parser.parse_args()

Folder=args.MosedepthOutPutFolder

# In[80]:

files = glob.glob(Folder+"*thresholds.bed.gz")


for i in files :
    
    File=pd.read_csv(i,sep="\t")
    
    File["TotalBases"]=File["end"]-File["start"]
    
    DropCol=['#chrom', 'start', 'end', 'region','TotalBases']
    Columns=File.columns
    Columns.drop(DropCol)
    
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].apply(pd.to_numeric)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].div(File.TotalBases, axis=0)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].mul(100)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].astype(int)
    
    File.to_csv(i[:-7]+"_percentage_bed.tsv",sep="\t",float_format='%.2f',index=False)


##This script will calculate mean depth at different threshold across samples
coverage_dict={}

files=glob.glob(Folder+"*_recal.mosdepth.region.dist.txt")

for file in files:
    df =pd.read_csv(file,sep="\t",header=None)
    df=df[ (df[0]=="total") & (df[1]<101)].sort_values(by=1).iloc[1:,]
    sample=file.split("/")[-1].replace("_recal.mosdepth.region.dist.txt","")
    coverage_dict[sample]=df[2].to_list()

average_depth=pd.DataFrame(coverage_dict).T
average_depth.columns=list(range(1, 101))
average_depth.columns=[str(x)+"_X" for x in average_depth ]
average_depth=average_depth.reset_index()
average_depth=average_depth.rename(columns={"index":"sample"})
average_depth=average_depth.sort_values(by="20_X")


average_depth.to_csv(Folder+"Average_Exome_coverage_at_differentThreshold.csv",index=None)
