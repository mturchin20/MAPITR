#!/bin/sh

#20200624

mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations 
mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR 

#zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz

#scp -p mturchin@desktop4.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/.
#gunzip /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz

cat /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz 

#R -q -e "set.seed(973459); Data1 <- read.table(\"/users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
R -q -e "set.seed(973459); Data1 <- read.table(\"/users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); \
Y <- runif(nrow(Data1)); \ 
Pathways <- c(); \
Pathways.Full <- sample(1:ncol(Data1), 6000); \
Count1 <- 1; for (i in 1:30) { \
	Pathways <- rbind(Pathways, Pathways.Full[Count1:(Count1+29)]); \
	Count1 <- Count1 + 29; \
}; \
SNPs.Pathways <- c(); \
SNPs.Genome <- c(); \
Genome.SNPs <- 1:ncol(Data1); \
for (j in 1:4) { \
	SNPs.Pathways <- rbind(SNPs.Pathways, sample(Pathways[1,], 10)); \
	SNPs.Genome <- rbind(SNPs.Genome, sample(Genome.SNPs[! Genome.SNPs %in% Pathways[1,]], 10)); \
}; \
for (k in 1:4) { \
	for (l in 1:10) { \
		Y <- Y + .5 * (SNPs.Pathways[k,l] * SNPs.Genome[k,l]); \
	}; \
}; \
write.table(Y, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \
write.table(Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
write.table(SNPs.Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
write.table(SNPs.Genome, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
"

#30 pathways
#4 true
#each have 200 SNPs 
#10 real SNPs each interacting with 10 other SNPs in the genome


