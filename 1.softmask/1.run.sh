#softmask
nohup RepeatMasker -species insects -xsmall  -parallel 4 -dir pacbioSim_mask_out -gff ../0.genome/Drosophila_simulans_Genome.fasta > pacbioSim_mask.log 2>&1 &

#2bit file
gzip pacbioSim_mask_out/Drosophila_simulans_Genome.fasta.masked 

faToTwoBit pacbioSim_mask_out/Drosophila_simulans_Genome.fasta.masked.gz pacbioSim_mask_out/Drosophila_simulans_Genome.2bit 

twoBitInfo pacbioSim_mask_out/Drosophila_simulans_Genome.2bit stdout|sort -k2,2nr >pacbioSim_mask_out/Drosophila_simulans_Genome.chrom.sizes
