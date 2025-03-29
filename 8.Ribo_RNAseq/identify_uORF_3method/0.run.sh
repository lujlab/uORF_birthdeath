#1.ribotish
for i in ../bam/*sortedByCoord.psite.tag.bam; do samtools index $i ;done
for i in ../bam/*sortedByCoord.psite.tag.bam; do j=`basename $i Aligned.sortedByCoord.psite.tag.bam`; echo "python ~/lcl/software/psite_ZH/src/ribotish_predict.py -p4 --framebest -b $i -g /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf -f /gpfs2/liucl/lcl/dataset/dm6/dm6.04/dmel-all-chromosome-r6.04.fasta -o ribotish/${j}.psite.ribotish_pred.txt --allresult ribotish/${j}.psite.ribotish_allresult.txt" >> 1.ribotish.sh; done

#2.ribocode
for i in ../bam/*toTranscriptome.psite.tag.bam; do j=`basename $i Aligned.toTranscriptome.psite.tag.bam`; echo "python ~/lcl/software/psite_ZH/src/ribocode_config.py $i ${i%%Aligned*}.psite.log > ribocode_config/${j}.pre_config" >> 2.1.ribocode_config.sh; done

#merge replicates
cat dmel_male_body1_ribo.pre_config dmel_male_body2_ribo.pre_config >dmel_male_body_ribo.pre_config
cat dmel_female_body1_ribo.pre_config dmel_female_body2_ribo.pre_config >dmel_female_body_ribo.pre_config
cat Dunn_em02_rep1_ribo.pre_config  Dunn_em02_rep2_ribo.pre_config > Dunn_em02_ribo.pre_config
for i in Samuels_GSC_AHS*_rep1_ribo.pre_config;do j=`basename $i _rep1_ribo.pre_config`; cat $i ${j}_rep2_ribo.pre_config ${j}_rep3_ribo.pre_config ${j}_rep4_ribo.pre_config > ${j}_ribo.pre_config;done


for i in ribocode_config/*.pre_config; do j=`basename $i .pre_config`; echo "python ~/lcl/software/psite_ZH/src/RiboCode_psite.py -a /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/RiboCode_annot_BDGP6.32.52_main -c ribocode_config/${j}.pre_config -l no -g -o ribocode/${j}.psite" >> 2.ribocode.sh; done

#3.RibORF2

#get Asite sam
for i in ../bam/*sortedByCoord.psite.tag.bam;do j=`basename $i Aligned.sortedByCoord.psite.tag.bam`; echo "python  ~/lcl/software/psite_ZH/psite-main/psite/pbam_Asite.py -f sam -l 27 -u 34 /gpfs2/liucl/lcl/dataset/dm6/dm6.04/dmel-all-chromosome-r6.04.fasta $i ../bam/${j}.psite.gbt.pickle Asite/${j}_Asite.sam">>3.1.getAsite.sh;done
#merge replicates
for i in `ls Asite/Samuels_GSC_AHS*_rep1_ribo_Asite.sam`;do j=`basename $i _rep1_ribo_Asite.sam`; echo "samtools merge Asite/${j}_ribo_Asite.sam $i Asite/${j}_rep2_ribo_Asite.sam Asite/${j}_rep3_ribo_Asite.sam Asite/${j}_rep4_ribo_Asite.sam" >>3.2.merge.sh;done
echo "samtools merge Asite/dmel_female_body_ribo_Asite.sam Asite/dmel_female_body1_ribo_Asite.sam Asite/dmel_female_body2_ribo_Asite.sam" >>3.2.merge.sh
echo "samtools merge Asite/dmel_male_body_ribo_Asite.sam Asite/dmel_male_body1_ribo_Asite.sam Asite/dmel_male_body2_ribo_Asite.sam" >>3.2.merge.sh
echo "samtools merge Asite/Dunn_em02_ribo_Asite.sam Asite/Dunn_em02_rep1_ribo_Asite.sam Asite/Dunn_em02_rep2_ribo_Asite.sam" >>3.2.merge.sh


#nohup python  ~/lcl/software/psite_ZH/psite-main/psite/pbam_Asite.py -f sam -l 26 -u 34 /gpfs2/liucl/lcl/dataset/dm6/dm6.04/dmel-all-chromosome-r6.04.fasta ../psite/1.remap_new/bam/Kronja_mature_oocyte_ribo_Aligned.sortedByCoord.out.bam ../psite/1.remap_new/bam/Kronja_mature_oocyte_ribo.psite.gbt.pickle Kronja_mature_oocyte_ribo_Asite.sam &

#run ribORF.pl to predict ORFs

#perl ~/lcl/software/RibORF2.0/RibORF-master/RibORF.2.0/ribORF.pl -f Kronja_mature_oocyte_ribo_Asite.sam -c ORFannotate_out/candidateORF.genepred.txt -o out_ribORF
for i in Asite/*_ribo_Asite.sam;do j=`basename $i _ribo_Asite.sam`;echo "mkdir ribORF_old/${j}_riborf;perl ~/lcl/software/RibORF2.0/RibORF-master/RibORF.2.0/ribORF.pl -f $i -c ../../2.2.identify_orfs/RibORF2/ORFannotate_out/candidateORF.genepred.txt -o ribORF/${j}_riborf > ribORF/${j}_riborf/log.out 2>&1" >>3.3.riborf.sh ;done

# 4.ORF_classification
for i in ribocode/*.psite.txt; do j=`basename $i .psite.txt` ;echo "python ~/lcl/software/psite_ZH/src/ncorf_classifier3.py -k -m ribocode -p allorfs/${j}.ribocode $i /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf.info2 > allorfs/${j}.ribocode.log 2>&1 ">>4.1.class_ribocode.sh; done
for i in ribotish/*psite.ribotish_pred.txt; do j=`basename $i .psite.ribotish_pred.txt`; echo "python ~/lcl/software/psite_ZH/src/ncorf_classifier3.py -k -m ribotish -p allorfs/${j}.ribotish $i /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf.info2 >allorfs/${j}.ribotish.log 2>&1 ">>4.2.class_ribotish.sh; done

for i in `less samplelist`;do echo "python ~/lcl/software/psite_ZH/src/ncorf_classifier3.py -k -m riborf -p allorfs/${i}_ribo.riborf ribORF/${i}_riborf/repre.valid.pred.pvalue.parameters.txt /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.2.identify_orfs/psite/0.ann/Drosophila_melanogaster.BDGP6.32.52.main.gtf.info2 >allorfs/${i}.riborf.log 2>&1">>4.3.class_riborf.sh; done 


