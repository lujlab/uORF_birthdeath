mkdir scratch
mkdir scratch/tmp
chmod +x runLastz
gensub2 target.list query.list template jobList

nohup ParaFly -c jobList -CPU 10 &
csh chainJobs.csh
chainPreNet dm6.PacBioSim.all.chain.gz dm6.chrom.sizes PacBioSim.chrom.sizes dm6.PacBioSim.all.pre.chain
chainNet dm6.PacBioSim.all.pre.chain  dm6.chrom.sizes PacBioSim.chrom.sizes dm6.prenet PacBioSim.prenet

netSyntenic dm6.prenet dm6.net
netSyntenic PacBioSim.prenet PacBioSim.net

netToAxt dm6.net dm6.PacBioSim.all.pre.chain dm6.2bit PacBioSim.2bit dm6_PacBioSim.axt
axtSort dm6_PacBioSim.axt dm6_PacBioSim_sort.axt

axtToMaf dm6_PacBioSim_sort.axt dm6.chrom.sizes PacBioSim.chrom.sizes dm6_PacBioSim.maf -tPrefix=dm6. -qPrefix=PacBioSim.
