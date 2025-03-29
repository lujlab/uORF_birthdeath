#psite_adjusted
cp /lustre/user/lulab/luoshq/syq/reference/dmel-all-r6.04.clean.gtf_rois.txt .

#1_psite
for i in `ls ../mel/map2genome/*riboAligned.sortedByCoord.out.bam` ; do j=`basename $i Aligned.sortedByCoord.out.bam`; echo "samtools index $i ; psite --count_files $i --countfile_format BAM --min_length 27 --max_length 34 --require_upstream /lustre/user/lulab/luoshq/syq/reference/dmel-all-r6.04.clean.gtf_rois.txt plastid/${j}.plastid">>1_psite.sh;done

#2_Psite_cnt
for i in `ls ../mel/map2genome/*riboAligned.sortedByCoord.out.bam` ; do j=`basename $i Aligned.sortedByCoord.out.bam`; echo "python bam_coverage.py -l 27 -u 34 -v plastid/${j}.plastid_p_offsets.txt $i Psite/${j}_ ;/lustre/user/lulab/doushq/wuxk/software/bigwigtobedgraph/bigWigToBedGraph Psite/${j}__fw.bw Psite/${j}__fw.bedgraph ; /lustre/user/lulab/doushq/wuxk/software/bigwigtobedgraph/bigWigToBedGraph Psite/${j}__rc.bw Psite/${j}__rc.bedgraph" >> 2_Psite_cnt.sh;done
for i in `ls ../mel/map2genome/*riboAligned.sortedByCoord.out.bam` ; do j=`basename $i Aligned.sortedByCoord.out.bam`; cat Psite/${j}__fw.bedgraph|awk -vOFS='\t' '{print $1,$2,$3,".",$4,"+"}' > Psite/${j}_riboseq_strand.bedgraph; cat Psite/${j}__rc.bedgraph|awk -vOFS='\t' '{print $1,$2,$3,".",$4,"-"}' >> Psite/${j}_riboseq_strand.bedgraph;done
for i in `ls Psite/*_riboseq_strand.bedgraph`; do j=`basename $i _riboseq_strand.bedgraph`; echo "sort -T ~/tmp -k1,1 -k2,2n $i > Psite/${j}_sorted_riboseq_strand.bedgraph" >>2_2_sort.sh;done

#librarysize
for i in `ls Psite/*_sorted_riboseq_strand.bedgraph`; do j=`basename $i _ribo_sorted_riboseq_strand.bedgraph`;awk 'BEGIN{sum=0}{sum += $5}END{print sum}' $i > lib_size/${j}_library_size;done

#3.collapse_exon_psite.sh
#sed 's/chr//g' ../../../8.new_27_mfa/PacBioSim/ann/dmel-all-r6.04.CDS.bed >dmel-all-r6.04.CDS.bed

mkdir ribo_pnt_coverage
for i in `ls Psite/*_sorted_riboseq_strand.bedgraph`; do j=`basename $i _ribo_sorted_riboseq_strand.bedgraph`;echo "bedtools map -a ../mel/r6.04.transcript.exon.win1nt.bed -b $i -s -c 5 -o sum -sorted -g ../mel/dmel-all-chromosome-r6.04.fasta.genomefile | sort -T ~/tmp -k4,4 -k5,5n | python ../bed_groupby.py > ribo_pnt_coverage/${j}.Psite.exon.bed" >>3_collapse_exon_psite.sh;done

#bedtools map -a ../mel/r6.04.transcript.exon.win1nt.bed -b Psite/dmel_em_0_2h_ribo_sorted_riboseq_strand.bedgraph -s -c 5 -o sum -sorted -g ../mel/dmel-all-chromosome-r6.04.fasta.genomefile|less

#4_contate_exon_psite.sh
# collapse exon cov to transcript coverage
for i in `ls ribo_pnt_coverage/*exon.bed`; do echo "python ../concatenate_exon.py $i ../mel/r6.04.transcript.exon.ordered.bed > ${i%%.exon.bed}.trans.bed">> 4_contate_exon_psite.sh;done

#5_sum_cds_reads.sh
mkdir all_cds_cnt
awk -vOFS='\t' '$2=="CDS" {print $1,$3,$4}' /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/0.ann/tr-r6.04-cds_ut_tr-pos.txt >all_cds_cnt/all_cds_tpos.txt
#python sum_pnt_cov.py mrna_pnt_coverage/mrna01_elife_repA.weighted.trans.bed all_cds_cnt/all_cds_tpos.txt > all_cds_cnt/mrna01_elife_repA_cds_cnt.bed
for i in `ls ribo_pnt_coverage/*.trans.bed`; do j=`basename $i .trans.bed`; echo "python ../sum_pnt_cov.py $i all_cds_cnt/all_cds_tpos.txt >all_cds_cnt/${j}_cds_cnt.bed">>5_sum_cds_reads.sh;done

#6_sum_uorf_reads.sh
mkdir all_uorf_cnt/
awk -vOFS='\t' '{print $1,$2+1,$3}' ../find_nonoverlaped_uORF_region/dmel_uORF_rmCDS_tr.bed > all_uorf_cnt/all_uorf_tpos.txt
#python sum_pnt_cov.py mrna_pnt_coverage/mrna01_elife_repA.weighted.trans.bed all_orf_cnt/all_orfs_tpos.txt > all_orf_cnt/mrna01_elife_repA_orf_cnt.bed
for i in `ls ribo_pnt_coverage/*.trans.bed`; do j=`basename $i .trans.bed`; echo "python ../sum_pnt_cov.py $i all_uorf_cnt/all_uorf_tpos.txt >all_uorf_cnt/${j}_uorf_cnt.bed">>6_sum_uorf_reads.sh;done

