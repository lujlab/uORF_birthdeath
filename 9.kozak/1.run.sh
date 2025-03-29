#sim PacBioSim
# kozak sequence
cat ../2.sim_coor/dsim_uORF_nogap.bed.corrected|awk 'OFS="\t"{if($2>=6 && $3<$7){print $1,$2-6,$2+4}}'|bedtools getfasta -fi ../2.sim_coor/dmel-all-r6.04.gtf.UTR5.12col.chr.newSim.nogap.fasta -bed - -fo -|perl -wpe 's/>(\S+):([0-9]+)\-([0-9]+)\n/$1\t$2\t$3\t/g'|awk 'OFS="\t"{print $1,$2+6,$4}'|perl add.intersect.2.pl - ../2.sim_coor/dsim_uORF_nogap.bed.corrected 1 2 3 1 2 N > dsim_uORF_nogap.bed.corrected.addkozak
#kozak score
cat dsim_uORF_nogap.bed.corrected.addkozak|awk '{if($11!="N"){print $11}}'|Rscript kozak.R /gpfs/user/liuchx/find_uORF/dmel6.cds_kozak.txall.txt -|perl ../2.sim_coor/add.intersect.pl - dsim_uORF_nogap.bed.corrected.addkozak 1 2 11 NA > dsim_uORF_nogap.bed.corrected.addkozakscore
cat dsim_uORF_nogap.bed.corrected.addkozak|awk '{if($11!="N"){print $11}}'|Rscript kozak.R /gpfs/user/liuchx/find_uORF/dsim1.cds_kozak.txall.txt -|perl ../2.sim_coor/add.intersect.pl - dsim_uORF_nogap.bed.corrected.addkozak 1 2 11 NA > dsim_uORF_nogap.bed.corrected.addkozakscore2



#dm6
less ../../8.new_27_mfa/PacBioSim/dm6_PacBioSim_27species_5UTR.fa|awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }'|grep 'dm6'|sed -e 's/-//g' -e 's/\t/\n/g' > dmel-all-r6.04.UTR5.12.maf.dm6.nogap.fa
less ../../8.new_27_mfa/PacBioSim/dm6_PacBioSim_27species_5UTR.fa|awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }'|grep 'dm6'|sed -e 's/\t/\n/g' > dmel-all-r6.04.UTR5.12.maf.dm6.gap.fa

#生成 dm6 uORF的 bed
less ../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update.csv|tail -n +2|sed 's/,/\t/g'|cut -f 1-3|awk '$2+$3!=0'|perl -wpe 's/\(\d+\).+//g'|sed 's/_/\t/g'|awk -v OFS="\t" '{print "dm6."$1,0,$2+1,$1"_"$2,"*","+"}' > dm6_uA_gap.bed
bedtools getfasta -s -fi dmel-all-r6.04.UTR5.12.maf.dm6.gap.fa -bed dm6_uA_gap.bed -fo test_seq 

less test_seq |awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }'|awk -v OFS="\t" '{gsub(/-/,"",$2); print $1,length($2)}'|perl -wpe 's/>(dm6.FBtr\d+):0-(\d+)\(.\)(.+)/$1\t$2\t$3/g'|awk -v OFS="\t" '{if($3>0) print $1,$3-1,$3+2,$1"_"$2"_"$3,"*","+"}' > dm6_uATG_nogap.bed

# upper case
cat dmel-all-r6.04.UTR5.12.maf.dm6.nogap.fa|sed 's/[a-z]/\u&/g'|perl -wpe 's/DM6.FBTR/dm6.FBtr/g' > dmel-all-r6.04.gtf.UTR5.12col.chr.dm6.nogap.fasta
# correct ungapped position with -1 0
cat dm6_uATG_nogap.bed|perl -wpe 's/_([0-9]+)_([0-9]+)\t\*/\t$1\t$2/g'|awk 'OFS="\t"{if($2<0){print $1,0,3,$4"_"$5,1,$7} else{print $1,$2,$3,$4"_"$5,$6,$7}}'|perl -wpe 's/\tdm6\./\t/g' > dm6_uATG_nogap.bed.tmp
# correct out of range
cat dmel-all-r6.04.gtf.UTR5.12col.chr.dm6.nogap.fasta|perl -wpe 's/>(\S+)\n/$1\t/g'|awk 'OFS="\t"{print $1,length($2)}'|perl ../2.sim_coor/add.intersect.pl - dm6_uATG_nogap.bed.tmp 1 2 1 0|awk 'OFS="\t"{if($3>$7){$3=$7;print $0} else{print}}'|perl ../2.sim_coor/add.intersect.pl ../2.sim_coor/uATG_matrix_dm6_FBtr_pos.bed - 7 4 4 NA|perl ../2.sim_coor/add.intersect.pl ../2.sim_coor/uATG_matrix_dm6_FBtr_pos.bed - 7 5 4 NA > dm6_uATG_nogap.bed.corrected
# check uATG
cat dm6_uATG_nogap.bed.corrected|grep "dm6:1"|bedtools getfasta -fi dmel-all-r6.04.gtf.UTR5.12col.chr.dm6.nogap.fasta -bed - -fo -|grep -v ">"|sort|uniq -c
# get seq from uATG & extend
cat dm6_uATG_nogap.bed.corrected|awk 'OFS="\t"{print $1,$2,$7}'|bedtools getfasta -fi dmel-all-r6.04.gtf.UTR5.12col.chr.dm6.nogap.fasta -bed - -fo -|perl -wpe 's/>(\S+):([0-9]+)\-([0-9]+)\n/$1\t$2\t$3\t/g'|perl -wpe 's/([ACGTN]{3})/$1,/g'|perl -wpe 's/,(TAA|TGA|TAG)\S+/,$1/g'|perl -wpe 's/,\S{1,2}$//g'|perl -wpe 's/,//g'|awk 'OFS="\t"{print $1,$2,$2+length($4),$4,$3}'|paste - dm6_uATG_nogap.bed.corrected|awk 'OFS="\t"{print $1,$2,$3,$9,$10,$11,$12,$13,$14,$4}' > dm6_uORF_nogap.bed.corrected

# kozak sequence
cat dm6_uORF_nogap.bed.corrected|awk 'OFS="\t"{if($2>=6 && $3<$7){print $1,$2-6,$2+4}}'|bedtools getfasta -fi dmel-all-r6.04.gtf.UTR5.12col.chr.dm6.nogap.fasta -bed - -fo -|perl -wpe 's/>(\S+):([0-9]+)\-([0-9]+)\n/$1\t$2\t$3\t/g'|awk 'OFS="\t"{print $1,$2+6,$4}'|perl add.intersect.2.pl - dm6_uORF_nogap.bed.corrected 1 2 3 1 2 N > dm6_uORF_nogap.bed.corrected.addkozak
#kozak score
cat dm6_uORF_nogap.bed.corrected.addkozak|awk '{if($11!="N"){print $11}}'|Rscript kozak.R /gpfs/user/liuchx/find_uORF/dmel6.cds_kozak.txall.txt -|perl ../2.sim_coor/add.intersect.pl - dm6_uORF_nogap.bed.corrected.addkozak 1 2 11 NA > dm6_uORF_nogap.bed.corrected.addkozakscore




####
less ../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update.csv|tail -n +2|sed 's/,/\t/g'|cut -f 1-3|awk '$2+$3!=0'>uORF_id_tmp
paste uORF_id_tmp dm6_uORF_nogap.bed.corrected.addkozakscore >dm6_addkozakscore_tmp
#merge kozakscore files of dm6 and dsim, merge_kozakscore_files.R, ->uORF_kozakscore.txt 
less ../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update.csv|tail -n +2|sed 's/,/\t/g'|cut -f 1-3,30-31,33|awk '$2+$3!=0'>tmp1
less ../../8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update.csv|tail -n +2|sed 's/,/\t/g'|cut -f 1-3|awk '$2+$3!=0'|perl -wpe 's/\(\d+\)//g' |sed 's/_/\t/g' |awk '{print $1"_"$2+1}' |paste - tmp1 >uORF_id_tmp2


cut -f 4,5 dm6_uORF_nogap.bed.corrected|sed 's/_/\t/g' |awk -v OFS='\t' '{print $1"_"$2,$1"_"$3}'>dm6_uORF_id
cut -f 4,5 dsim_uORF_nogap.bed.corrected.addkozak|sed 's/_/\t/g' |awk -v OFS='\t' '{print $1"_"$2,$1"_"$3}'>dsim_PacBioSim_uORF_id

