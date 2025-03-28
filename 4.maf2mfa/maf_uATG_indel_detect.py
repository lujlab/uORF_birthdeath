#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
#by LIUchenlu 2023/2/13
#usage maf_uATG_indel_detect.py  input_uATG_coordinate input_maf_file temp_prefix outfile
###python cheat sheet
import sys
import os
#import numpy as np
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO


def load_data(file_path):
    f=open(file_path,"r") 
    data=[]
    for line in f.readlines():
        lines=line.replace(" \n","").replace("\n","").split("\t")
        data.append(lines)
    f.close()
    return data


def ATG_maf_gap_detect(trid,pos,species,species_list,seq_record):
    p1=pos
    tmp_id=species+"."+trid
    d=dict(zip(species_list,["NA"]*len(species_list)))
    for seq in seq_record: #find gaps in ATG in this species
        if seq.id==tmp_id:
            sequence=str(seq.seq).replace("a","A").replace("t","T").replace("c","C").replace("g","G")
            ATG_p1=sequence[p1:len(sequence)].index("A")+p1
            #ATG_p2=sequence[p1:len(sequence)].index("T")+p1
            ATG_p3=sequence[p1:len(sequence)].index("G")+p1
            gap_1=sequence[ATG_p1:ATG_p3+1].count("-") #number of gaps in ATG
            #gap_1_extended=extended_gap(seq=sequence,ATG_p1=ATG_p1,ATG_p3=ATG_p3,gap=gap_1)
            d[species]="0"
    for seq in seq_record:
        if seq.id!=tmp_id:
            s=seq.id.replace(trid,"").replace(".","")
            gap_2=seq.seq[ATG_p1:ATG_p3+1].count("-")
            if gap_2==gap_1:
                d[s]="0"
            else : #identify gap length
                relative_gap_len=extended_gap(seq1=sequence,seq2=seq.seq,ATG_p1=ATG_p1,ATG_p3=ATG_p3)
                d[s]=str(relative_gap_len)
    return d

def extended_gap(seq1,seq2,ATG_p1,ATG_p3):
    #5' 
    seq1_5=''.join(reversed(seq1[0:ATG_p1])).replace("A","N").replace("T","N").replace("C","N").replace("G","N").replace("a","N").replace("t","N").replace("c","N").replace("g","N").replace("n","N")
    seq2_5=''.join(reversed(seq2[0:ATG_p1])).replace("A","N").replace("T","N").replace("C","N").replace("G","N").replace("a","N").replace("t","N").replace("c","N").replace("g","N").replace("n","N")
    seq1_5=list(seq1_5.replace("N","1").replace("-","0"))
    seq2_5=list(seq2_5.replace("N","1").replace("-","0"))
    seq1_5 = [ int (i) for i in seq1_5 ]
    seq2_5 = [ int (i) for i in seq2_5 ]
    seq12_5 = []
    for x, y in zip(seq1_5, seq2_5):
        seq12_5.append(x + y)
    seq12_5.append(2)
    r1=seq12_5.index(2)

    #3'
    seq1_3=str(seq1[ATG_p3+1:len(seq1)]).replace("A","N").replace("T","N").replace("C","N").replace("G","N").replace("a","N").replace("t","N").replace("c","N").replace("g","N").replace("n","N")
    seq2_3=str(seq2[ATG_p3+1:len(seq2)]).replace("A","N").replace("T","N").replace("C","N").replace("G","N").replace("a","N").replace("t","N").replace("c","N").replace("g","N").replace("n","N")
    seq1_3=list(seq1_3.replace("N","1").replace("-","0"))
    seq2_3=list(seq2_3.replace("N","1").replace("-","0"))
    seq1_3 = [ int (i) for i in seq1_3 ]
    seq2_3 = [ int (i) for i in seq2_3 ]
    seq12_3 = []
    for x, y in zip(seq1_3, seq2_3):
        seq12_3.append(x + y)
    seq12_3.append(2)
    r2=seq12_3.index(2)

    gap=seq1[ATG_p1-r1:ATG_p3+r2+1].count("-")-seq2[ATG_p1-r1:ATG_p3+r2+1].count("-")
    return gap


if __name__ == '__main__':
    uATG_coordinate = sys.argv[1]  #uATG_coordinate="/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/uORF_matrix_triCas2_update.tsv"
    fasta_align = sys.argv[2] #fasta_align = "/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dmel-all-r6.04.gtf.UTR5.12col.chr.maf"
    tmp_prefix = sys.argv[3] #tmp_prefix="tmp1"
    outfile = sys.argv[4]

    f1=load_data(uATG_coordinate)
    head=f1[0]
    species_list=f1[0][2:30]
    f1=pd.DataFrame(f1,columns=f1[0])
    f1.drop([0],inplace=True)
    tr_list=list(set(list(f1["trid"])))

    #print head
    print("\t".join(head)+"\t"+"_gap\t".join(species_list)+"_gap", file=open(outfile, "a"))

    for i in tr_list:
        cmd="grep "+i+" "+fasta_align+" |sed 's/>//g' |seqkit grep -f - "+fasta_align+">"+tmp_prefix+i+".fa"
        os.system(cmd)
        seq_record=list(SeqIO.parse(tmp_prefix+i+".fa", "fasta"))
        tr_uATG_list=f1[ f1.trid == i].reset_index(drop=True)
        for j in range(0,len(tr_uATG_list["maf_pos"])):
            #find position
            p=tr_uATG_list["maf_pos"][j].replace("("," ").replace(")","").split(" ")
            p1=int(p[0])
            #find the first species with ATG
            x=tr_uATG_list.loc[j].tolist()[2:30]
            species_index=x.index("1")
            s=species_list[species_index]
            tmp_id=s+"."+i
            #flag
            #flag_indel="NA" 
            #find gaps in ATG in this species
            d=ATG_maf_gap_detect(trid=i,pos=p1,species=s,species_list=species_list,seq_record=seq_record)
            out1="\t".join(list(tr_uATG_list.loc[j]))
            out2="\t".join(list(d.values()))

            print(out1+"\t"+out2, file=open(outfile, "a"))
        os.system("rm "+tmp_prefix+i+".fa")
