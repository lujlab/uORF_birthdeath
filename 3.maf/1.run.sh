multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chr2L.maf ../chr2L.maf 1 chr2L_out1.maf chr2L_out2.maf all> dm6_PacBioSim_27species.chr2L.maf #the ../chr2L.maf file is downloaded from ucsc 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chr2R.maf ../chr2R.maf 1 chr2R_out1.maf chr2R_out2.maf all> dm6_PacBioSim_27species.chr2R.maf 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chr3L.maf ../chr3L.maf 1 chr3L_out1.maf chr3L_out2.maf all> dm6_PacBioSim_27species.chr3L.maf 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chr3R.maf ../chr3R.maf 1 chr3R_out1.maf chr3R_out2.maf all> dm6_PacBioSim_27species.chr3R.maf 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chr4.maf ../chr4.maf 1 chr4_out1.maf chr4_out2.maf all> dm6_PacBioSim_27species.chr4.maf 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chrX.maf ../chrX.maf 1 chrX_out1.maf chrX_out2.maf all> dm6_PacBioSim_27species.chrX.maf 
multiz ../../4.maffilter/pacbioSim/dm6_PacBioSim.chrY.maf ../chrY.maf 1 chrY_out1.maf chrY_out2.maf all> dm6_PacBioSim_27species.chrY.maf 
maffilter input.file=chr2L_out2.maf input.file.compression=none output.log=2L.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr2L)",Output"(file=chr2L_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chr2R_out2.maf input.file.compression=none output.log=2R.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr2R)",Output"(file=chr2R_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chr3L_out2.maf input.file.compression=none output.log=3L.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr3L)",Output"(file=chr3L_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chr3R_out2.maf input.file.compression=none output.log=3R.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr3R)",Output"(file=chr3R_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chr4_out2.maf input.file.compression=none output.log=4.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr4)",Output"(file=chr4_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chrX_out2.maf input.file.compression=none output.log=X.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chrX)",Output"(file=chrX_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chrY_out2.maf input.file.compression=none output.log=Y.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chrY)",Output"(file=chrY_out_add.maf, compression=none, mask=yes)"
maffilter input.file=chr2L_out1.maf input.file.compression=none output.log=2L.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr2L)",Output"(file=chr2L_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chr2R_out1.maf input.file.compression=none output.log=2R.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr2R)",Output"(file=chr2R_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chr3L_out1.maf input.file.compression=none output.log=3L.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr3L)",Output"(file=chr3L_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chr3R_out1.maf input.file.compression=none output.log=3R.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr3R)",Output"(file=chr3R_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chr4_out1.maf input.file.compression=none output.log=4.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chr4)",Output"(file=chr4_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chrX_out1.maf input.file.compression=none output.log=X.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chrX)",Output"(file=chrX_out1_add.maf, compression=none, mask=yes)"
maffilter input.file=chrY_out1.maf input.file.compression=none output.log=Y.log maf.filter=SelectChr"(ref_species=dm6,chromosome=chrY)",Output"(file=chrY_out1_add.maf, compression=none, mask=yes)"
cat dm6_PacBioSim_27species.chr2L.maf chr2L_out1_add.maf chr2L_out_add.maf >dm6_PacBioSim_all.chr2L.maf
cat dm6_PacBioSim_27species.chr2R.maf chr2R_out1_add.maf chr2R_out_add.maf >dm6_PacBioSim_all.chr2R.maf
cat dm6_PacBioSim_27species.chr3L.maf chr3L_out1_add.maf chr3L_out_add.maf >dm6_PacBioSim_all.chr3L.maf
cat dm6_PacBioSim_27species.chr3R.maf chr3R_out1_add.maf chr3R_out_add.maf >dm6_PacBioSim_all.chr3R.maf
cat dm6_PacBioSim_27species.chr4.maf chr4_out1_add.maf chr4_out_add.maf >dm6_PacBioSim_all.chr4.maf
cat dm6_PacBioSim_27species.chrX.maf chrX_out1_add.maf chrX_out_add.maf >dm6_PacBioSim_all.chrX.maf
cat dm6_PacBioSim_27species.chrY.maf chrY_out1_add.maf chrY_out_add.maf >dm6_PacBioSim_all.chrY.maf
