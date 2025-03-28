sed 's/-//g' ../dm6_PacBioSim_27species_5UTR.fa >dm6_PacBioSim_27species_5UTR_nongap.fa
nohup python getfasta_PacBioSim_27species_merged.py &
sed 's/-//g' dm6_PacBioSim_27species_merged_5UTR.fa >dm6_PacBioSim_27species_merged_5UTR_nongap.fa
 ~/software/ucsc-kent/faCount -dinuc dm6_PacBioSim_27species_merged_5UTR_nongap.fa > dm6_PacBioSim_27species_merged_5UTR_nongap.dinuc.txt




