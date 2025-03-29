library(data.table)
library(UpSetR)


#uORF matrix and uniq position
uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix$geno_pos=paste0(uORF_matrix$seqname,"_",uORF_matrix$gene_pos)
uORF_matrix$uorf_id=paste0(uORF_matrix$transcriptID,"_",as.character(uORF_matrix$position_withoutgap+1))
uORF_matrix_unique=uORF_matrix[!duplicated(uORF_matrix$geno_pos),]
uORF_matrix_unique=uORF_matrix_unique[dm6==1,]
uATG_position=uORF_matrix_unique$geno_pos


#gene id and transcript id
tr_gn=fread("/data/r6.04.gene_tr",header=F)
colnames(tr_gn)=c("gn","tx_name")


sample <- fread('less samplelist21', header = F)
for (i in sample$V1){
  tmp_riborf=fread(paste0("allorfs/",i,"_ribo.riborf_processed.tsv"),header=T,sep='\t')
  tmp_riborf=tmp_riborf[,c("tx_name","start_codon","orf_type","tstart")]
  tmp_riborf$method="riborf"
  tmp_ribocode=fread(paste0("allorfs/",i,"_ribo.ribocode_processed.tsv"),header=T,sep='\t')
  tmp_ribocode=tmp_ribocode[,c("tx_name","start_codon","orf_type","tstart")]
  tmp_ribocode$method="ribocode"
  tmp_ribotish=fread(paste0("allorfs/",i,"_ribo.ribotish_processed.tsv"),header=T,sep='\t')
  tmp_ribotish=tmp_ribotish[,c("tx_name","start_codon","orf_type","tstart")]
  tmp_ribotish$method="ribotish"
  tmp=rbind(tmp_riborf,tmp_ribocode,tmp_ribotish)
  tmp=tmp[start_codon=="ATG",]
  tmp$uorf_id=paste0(tmp$tx_name,"_",tmp$tstart)
  assign(paste0(i,"_result"),tmp)
}

###stat uORF number (unique uATG pos)
uORF_stat=matrix(ncol=21,nrow=4)
colnames(uORF_stat) <- sample$V1
rownames(uORF_stat) <- c("riborf","ribocode","ribotish","combined")
upset_data_uORF=matrix(ncol=21,nrow=4)
upset_data_uORF=data.table(matrix(ncol = 21, nrow = 36565))
colnames(upset_data_uORF)=sample$V1
upset_data_uORF$geno_pos=uATG_position
upset_data_uORF[is.na(upset_data_uORF)]=0

total_riborf=data.table()
total_ribocode=data.table()
total_ribotish=data.table()
for (i in sample$V1){
  tmp=get(paste0(i,"_result"))
  #tmp_uorf=tmp[orf_type=="uORF" |orf_type=="N_extension" |orf_type=="uoORF" ,c("uorf_id")]
  uorf_id=tmp[,c("uorf_id","method")]
  uorf_id_pos=merge(uorf_id,uORF_matrix,by="uorf_id")
  #number of riborf identified uORF
  uORF_stat["riborf",i]=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="riborf",]$geno_pos,.N]
  uORF_stat["ribocode",i]=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="ribocode",]$geno_pos,.N]
  uORF_stat["ribotish",i]=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="ribotish",]$geno_pos,.N]
  uORF_stat["combined",i]=uORF_matrix_unique[geno_pos %in% uorf_id_pos$geno_pos,.N]
  upset_data_uORF[geno_pos %in% uorf_id_pos$geno_pos,i]=1
  total_riborf_tmp=data.table(geno_pos=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="riborf",]$geno_pos,]$geno_pos,sample=i)
  total_riborf=rbind(total_riborf,total_riborf_tmp)
  total_ribocode_tmp=data.table(geno_pos=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="ribocode",]$geno_pos,]$geno_pos,sample=i)
  total_ribocode=rbind(total_ribocode,total_ribocode_tmp)
  total_ribotish_tmp=data.table(geno_pos=uORF_matrix_unique[geno_pos %in% uorf_id_pos[method=="ribotish",]$geno_pos,]$geno_pos,sample=i)
  total_ribotish=rbind(total_ribotish,total_ribotish_tmp)
}

uORF_stat=data.frame(uORF_stat)
uORF_stat$method=rownames(uORF_stat)
fwrite(uORF_stat,"uORF_number.txt",sep='\t') #TableS4
###stat CDS number (unique gene id)
CDS_stat=matrix(ncol=21,nrow=4)
colnames(CDS_stat) <- sample$V1
rownames(CDS_stat) <- c("riborf","ribocode","ribotish","combined")
total_riborf=data.table()
total_ribocode=data.table()
total_ribotish=data.table()
for (i in sample$V1){
  tmp=get(paste0(i,"_result"))
  tmp=tmp[orf_type=="N_truncation" |orf_type=="CDS" ,c("tx_name","method")]
  tmp=merge(tmp,tr_gn,by="tx_name")
  #number of riborf identified uORF
  CDS_stat["riborf",i]=length(unique(tmp[method=="riborf",]$gn))
  CDS_stat["ribocode",i]=length(unique(tmp[method=="ribocode",]$gn))
  CDS_stat["ribotish",i]=length(unique(tmp[method=="ribotish",]$gn))
  CDS_stat["combined",i]=length(unique(tmp$gn))
  total_riborf_tmp=data.table(geno_pos=unique(tmp[method=="riborf",]$gn),sample=i)
  total_riborf=rbind(total_riborf,total_riborf_tmp)
  total_ribocode_tmp=data.table(geno_pos=unique(tmp[method=="ribocode",]$gn),sample=i)
  total_ribocode=rbind(total_ribocode,total_ribocode_tmp)
  total_ribotish_tmp=data.table(geno_pos=unique(tmp[method=="ribotish",]$gn),sample=i)
  total_ribotish=rbind(total_ribotish,total_ribotish_tmp)
}
CDS_stat=data.frame(CDS_stat)
CDS_stat$method=rownames(CDS_stat)
fwrite(CDS_stat,"CDS_number.txt",sep='\t') #TableS4



