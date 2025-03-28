import pandas as pd 
from Bio import SeqIO
import re

#species_list=["dm6","droSim1","droSec1","droYak3","droEre2","droBia2","droSuz1","droAna3","droBip2","droEug2","droEle2","droKik2","droTak2","droRho2","droFic2","droPse3","droPer1","droMir2","droWil2","droVir3","droMoj3","droAlb1","droGri2","musDom2","anoGam1","apiMel4","triCas2"]
species_list=["dm6","PacBioSim","droSim1","droSec1","droYak3","droEre2","droBia2","droSuz1","droAna3","droBip2","droEug2","droEle2","droKik2","droTak2","droRho2","droFic2","droPse3","droPer1","droMir2","droWil2","droVir3","droMoj3","droAlb1","droGri2","musDom2","anoGam1","apiMel4","triCas2"]
#species_list=["dm6","droSim1","droSec1","droYak3","droEre2","droBia2","droSuz1","droAna3","droBip2","droEug2","droEle2","droKik2","droTak2","droRho2","droFic2","droPse3","droPer1","droMir2","droWil2","droVir3","droMoj3","droAlb1","droGri2"]
records =list(SeqIO.parse("/data/dm6_PacBioSim_27species_5UTR.fa", "fasta"))

#select species in records
records=list(filter(lambda x: re.match(r"[\d\w]*",x.id).group() in species_list,records)) #there is difference of "filter" between python2 and python3.

#get mutialignment group by UTR area
import re
group_index_list=[]
UTR_id_list=[]
group_index_list.append(0)
index_num=0
while(index_num<len(records)-1):
    UTR_id=re.findall(r"FBtr.+",records[index_num].id)[0]
    while(index_num<len(records)-1 and UTR_id==re.findall(r"FBtr.+",records[index_num].id)[0]):
        index_num+=1
    UTR_id_list.append(UTR_id)
    if(index_num<len(records)-1):
        group_index_list.append(index_num)

uORF_matrix=pd.DataFrame(columns=species_list)
#find uATG in each mutialignment group
for i in range(0,len(group_index_list)):
    uATG_index_list={}         #store the position of uATG(start:end) in the mutialign
    uATG_genome_index_list={}  #store the position of uATG in genome without caculating gaps(muti-start:genome-start)
    if(i!=len(group_index_list)-1):
        temp=group_index_list[i+1]
    if(i==len(group_index_list)-1):
        temp=len(records)
    for j in range(group_index_list[i],temp):
        iterator=re.finditer(r"[Aa]-*[Tt]-*[Gg]",str(records[j].seq))
        for match in iterator:
            if(match.start() not in uATG_index_list):    #calculate the start position without any gap("-")
                uATG_genome_index_list[match.start()]=match.start()-records[group_index_list[i]].seq.count("-",0,match.start())
#if the new uORF found in species that are not dm6, and the start position of this uORF is with "-", we use next base's position as it.
                uATG_index_list[match.start()]=match.end()  #judge the situation that there are several gaps among uATG
    for uATG_index in uATG_index_list:
        binary_list=[]
        cursor1=0
        cursor2=group_index_list[i]
        if(i!=len(group_index_list)-1):
            cursor2_inf=group_index_list[i+1]-1
        if(i==len(group_index_list)-1):
            cursor2_inf=len(records)-1
        while(cursor1<len(species_list)):
            if(species_list[cursor1]==re.match(r"[\d\w]*",records[cursor2].id).group()):
                if(re.match(r"[Aa]-*[Tt]-*[Gg]",str(records[cursor2].seq)[uATG_index:uATG_index_list[uATG_index]])!=None):
                    binary_list.append("1")
                else:
                    binary_list.append("0")
                cursor1+=1
                if(cursor2<cursor2_inf):
                    cursor2+=1
            else:
                binary_list.append("0")
                cursor1+=1
        name=[]
        name.append(UTR_id_list[i]+"_"+str(uATG_index)+"("+str(uATG_genome_index_list[uATG_index])+")")
        newline=pd.DataFrame(dict(zip(species_list,binary_list)),index=name)
        uORF_matrix=uORF_matrix.append(newline,ignore_index=False)
uORF_matrix.to_csv("/results/uORF_matrix_dm6_PacBioSim_27species.csv")
