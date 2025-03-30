# uORF_brithdeath
I. Introduction
The purpose of these scripts is to investigate the birth, death, and evolutionary copensation of uORFs in Drosophila. Please refer to our article for specific experimental design and description.
More processed data and results are available at the Code Ocean Repository (10.24433/CO.0818531.v1).

II. Indication of script functions:
1.softmask: set repetitive sequences to lowercase (soft-masked) via RepeatMasker v4.1.1 (http://www.repeatmasker.org) for the PacBio D. simulans genome.

2.lastz: aligned the PacBio D. simulans genome to the D. melanogaster genome dm6 following UCSC guidelines.

3.maf: integrate this alignment into the multiple alignments of 27 species with Multiz.

4.maf2mfa: the 5' UTR sequence of each annotated transcript of D. melanogaster and its corresponding sequences in the other insects were extracted from the multiple alignments. All ATG triplets in the 5' UTRs of each species were identified.

5.GLOOM: infer the occurrence of gain and loss for each ATG based on maximum parsimony in multispecies analysis.

6.compensation: calculate the number of each class of uORFs and perform the bootstrap test for compensatory evolution of gain and loss events of uORFs.

7.mk_test: MK test for newly fixed uATG gain and loss events in D. melanogaster.

8.Ribo_RNAseq: Ribo-seq and mRNA-seq analysis.

9.kozak: calculate kozak scores for two species.

10.BLS: calculate BLS scores.

11.analysis_conserved_non-conserved_uORF：correlation between TE and conservation of uORFs

12.luciferase: experimental results and analysis of dual-luciferase report assays.

13.pop: uORF variants annotation and Tajima's D calculation.

III. Indication of figures in scripts:
Fig1: 6.compensation

Fig2: 11.analysis_conserved_nonconserve_new

Fig3: Fig3c-d: 6.compensation Fig3e: 7.mk_test

Fig4-5: 12.luciferase

Fig6: 13.pop

#supplementary figs: FigS1, S9-10, S25: 8.Ribo_RNAseq FigS2-8, S13-16: 5.GLOOM FigS12: 12.luciferase figS17-20: 11.analysis_conserved_nonconserve_new figS21: 6.compensation figS22: 7.mk_test figS23-24: 13.pop

V. License
Each file included in this capsule is licensed under the MIT License.
