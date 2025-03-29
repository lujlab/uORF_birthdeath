#get RNA seq coverage in Abd-B region

mkdir mrna_coverage
for i in `ls ../2.map3/mrna_weighted_bed/*.weighted.sorted.bed`; do j=`basename $i .weighted.sorted.bed`;echo "bedtools map -a Abd-B.1nt.bed -b $i -s  -c 5 -o sum -sorted -g /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel/dmel-all-chromosome-r6.04.fasta.genomefile >mrna_coverage/${j}_coverage">>mrna.sh;done

mkdir Psite_coverage

for i in `ls ../2.map3/coverage/*_sorted_riboseq_strand.bedgraph`;do j=`basename $i _sorted_riboseq_strand.bedgraph`;echo "bedtools map -a Abd-B.1nt.bed -b $i -s  -c 5 -o sum -sorted -g /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel/dmel-all-chromosome-r6.04.fasta.genomefile >Psite_coverage/${j}_coverage" >>Psite.sh;done

