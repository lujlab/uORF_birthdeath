import sys
import csv
import re
import os
import traceback
import pandas as pd
import tqdm as tqdm

gtf=pd.read_csv("../../../dmel-5UTR-r6.04.csv",index_col=0)
matrix=pd.read_csv("./uORF_matrix_dm6_PacBioSim_27species.csv",index_col=0)
matrix["transcriptID"]=matrix.index.map(lambda x: re.findall("FBtr[0-9]+",x)[0])
matrix["position_withoutgap"]=matrix.index.map(lambda x: re.findall("\((.*?)\)",x)[0])

def coordinate_trans(transcriptID,position):
    seqname="NA"
    geneID="NA"
    genepos=0
    genesymbol="NA"
    for index,row in gtf[gtf["transcript_id"]==transcriptID].iterrows():
        position=int(position)
        if( (position+1) > int(row["width"]) ):
            position-=int(row["width"])
            continue
        else:
            if (row["strand"] == "+"):
                genepos=position+row["start"]
            if (row["strand"] == "-"):
                genepos=row["end"]-position
            seqname=row["seqnames"]
            geneID=row["gene_id"]
            genesymbol=row["gene_symbol"]
            break
    return seqname,geneID,genepos,genesymbol
matrix["seqname"]=matrix.apply(lambda x:coordinate_trans(x["transcriptID"],x["position_withoutgap"])[0],axis=1)
matrix["geneID"]=matrix.apply(lambda x:coordinate_trans(x["transcriptID"],x["position_withoutgap"])[1],axis=1)
matrix["gene_pos"]=matrix.apply(lambda x:coordinate_trans(x["transcriptID"],x["position_withoutgap"])[2],axis=1)
matrix["genesymbol"]=matrix.apply(lambda x:coordinate_trans(x["transcriptID"],x["position_withoutgap"])[3],axis=1)
matrix.to_csv("/results/uORF_matrix_dm6_PacBioSim_27species_update.csv")
