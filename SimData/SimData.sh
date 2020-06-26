#!/bin/sh

#20200624

mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations 
mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR 

#zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz

#scp -p mturchin@desktop4.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/.
#gunzip /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz

#cat /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz 
cat /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit | head -n 5000 | perl -lane 'print join("\t", @F[0..99999]);' | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz 

#R -q -e "set.seed(973459); Data1 <- read.table(\"/users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
R -q -e "set.seed(973459); Data1 <- read.table(\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
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
Y.new <- Y; \
for (k in 1:4) { \
	for (l in 1:10) { \
		Y.new <- Y.new + .5 * (Data1[,SNPs.Pathways[k,l]] * Data1[,SNPs.Genome[k,l]]); \
	}; \
}; \
write.table(Y, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \
write.table(Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
write.table(SNPs.Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
write.table(SNPs.Genome, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
"

#R -q -e "
set.seed(973459); Data1 <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T); 
Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); 
Y <- rnorm(nrow(Data1)); 
Pathways <- c(); 
Pathways.Full <- sample(1:ncol(Data1), 6000); 
Count1 <- 1; for (i in 1:30) { 
	Pathways <- rbind(Pathways, Pathways.Full[Count1:(Count1+29)]); 
	Count1 <- Count1 + 29; 
}; 
SNPs.Pathways <- c(); 
SNPs.Genome <- c(); 
Genome.SNPs <- 1:ncol(Data1); 
for (j in 1:4) { 
	SNPs.Pathways <- rbind(SNPs.Pathways, sample(Pathways[1,], 10)); 
	SNPs.Genome <- rbind(SNPs.Genome, sample(Genome.SNPs[! Genome.SNPs %in% Pathways[1,]], 10)); 
}; 
Y.new <- Y; 
for (k in 1:4) { 
	for (l in 1:10) { 
		Y.new <- Y.new + .5 * (Data1[,SNPs.Pathways[k,l]] * Data1[,SNPs.Genome[k,l]]); 
	}; 
}; 
Pathways.Edits <- apply(Pathways, 1, function(x) { return(paste(x, collapse=",")); });
Pathways.Edits <- cbind(1:length(Pathways.Edits), Pathways.Edits);
Pathways.Edits <- cbind(rep("Pathway", nrow(Pathways.Edits)), Pathways.Edits);
Pathways.Edits.2 <- apply(Pathways.Edits[,1:2], 1, function(x) { return(paste(x, collapse="")); });
Pathways.Edits <- cbind(Pathways.Edits.2, Pathways.Edits[,3]); 
write.table(Y, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.Orig.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Y.new, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Pathways, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(Pathways.Edits, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(SNPs.Pathways, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(SNPs.Genome, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#"

#R -q -e"
library("devtools"); devtools::load_all();
set.seed(582724); 
X <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T);
Y <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", header=F);
Pathways <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", header=F);
SNPs.Pathways <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt", header=F);
SNPs.Genome <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt", header=F);

library("Rcpp"); library("RcppArmadillo"); library("RcppParallel"); library("doParallel");
sourceCpp("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/src/MAPITR.cpp");

library("devtools"); devtools::load_all();
set.seed(582724); 
X <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T);
Y <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.txt", header=F);
Pathways <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.Edits.txt", header=F);
SNPs.Pathways <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Pathways.txt", header=F);
SNPs.Genome <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Genome.txt", header=F);

library("Rcpp"); library("RcppArmadillo"); library("RcppParallel"); library("doParallel");
sourceCpp("/home/mturchin20/TempStuff/MAPITR/src/MAPITR.cpp")

Y2 <- cbind(Y,Y);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 



X.mean <- apply(X, 2, mean); X.sd <- apply(X, 2, sd); X.adj <- t((t(X)-X.mean)/X.sd); 
summary(lm(as.matrix(Y) ~ X.adj[,1252] + X.adj[,956] + X.adj[,1252] * X.adj[,956]))

##MAPITRmain <- function (Phenotype, Genotypes, Pathways, GRM_Grand = NULL, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) 
MAPITR.Results <- MAPITRmain(Y, X, Pathways);

#library("Rcpp"); library("RcppArmadillo"); 

MAPITRBase(as.matrix(Y2), as.matrix(t(X)), Pathways.Full[1:2], cores=1)
MAPITRBaseTest(as.matrix(Y), as.matrix(t(X)), Pathways.Full[1], cores=1)
MAPITRBaseTest(as.matrix(Y), Pathways.Full[1], cores=1)
Test2 <- 


#30 pathways
#4 true
#each have 200 SNPs 
#10 real SNPs each interacting with 10 other SNPs in the genome


