#### calculation branch length score ####
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/5.BLS/")
library(data.table)
library(ape) # https://cran.r-project.org/web/packages/ape/index.html

### these 3 BLS are not same even for a same uATG ###
tree_fly <- read.tree(file="/lustre/user/lulab/luoshq/syq/project/uorf_MK/gain_loss_syq/GLOOME/dm6.27way.commonNames.nh")
hh<-read.csv("../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",header = T,stringsAsFactors = F)
hh[,-4]->hh
names(hh)[3]<-"droSim1" #
list_fly <- apply(hh[,2:28],1,function(x){
  sp <- which(x>0)
  paste0(colnames(hh)[2:28][sp],collapse=",")
})
list_fly <- as.data.frame(list_fly)
colnames(list_fly) <- c("sp")
list_fly[,1] <- as.character(list_fly[,1])
# a long time for BLS calculation
list_fly$score <- 0
for(i in 1:nrow(list_fly)){
  list_fly$score[i] <- sum(keep.tip(tree_fly,unlist(strsplit(list_fly[i,1],",")))$edge.length)
}

total_score<-sum(tree_fly$edge.length) #11.399043
list_fly$BLS <- list_fly$score/total_score
names(list_fly)[1]<-"present_sp"
list_fly$present_sp<-gsub("droSim1","PacBioSim",list_fly$present_sp)

jj<-read.csv("../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",header = T,stringsAsFactors = F)
cbind(jj,list_fly)->hh2
write.table(hh2,"./uORF_matrix_triCas2_ATG_updata_BLS.txt",row.names = F,col.names = T,sep="\t",quote = F)
