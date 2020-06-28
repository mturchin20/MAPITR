#!/bin/sh

#20200624

mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations 
mkdir /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR 

#zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz

##scp -p mturchin@desktop4.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/.
#scp -p mturchin@desktop4.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz /home/mturchin20/TempStuff/MAPITR/SimData/. 
#gunzip /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz

#cat /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz 
cat /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit | head -n 5000 | perl -lane 'print join("\t", @F[0..99999]);' | R -q -e "set.seed(68361); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 10000); Rows <- sample(1:nrow(Data1), 2000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz 

##R -q -e "set.seed(973459); Data1 <- read.table(\"/users/mturchin/data/mturchin/InterPath/Analyses/Rnd2AdditiveMdls/Simulations/MAPITR/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
##R -q -e "set.seed(973459); Data1 <- read.table(\"/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz\", header=T); \
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
#write.table(Y, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \
#write.table(Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
#write.table(SNPs.Pathways, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
#write.table(SNPs.Genome, file=\"/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
"

#write.table(Y, file=\"/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \
#write.table(Pathways, file=\"/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
#write.table(SNPs.Pathways, file=\"/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Pathways.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 
#write.table(SNPs.Genome, file=\"/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Genome.txt\", quote=FALSE, row.names=FALSE, col.names=FALSE); \ 

#R -q -e "
#set.seed(973459); Data1 <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T); 
#R -q -e "
##set.seed(973459); Data1 <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz", header=T); 
##Data1 <- Data1[sample(1:nrow(Data1), 2000), sample(1:ncol(Data1), 10000)];
set.seed(973459); Data1 <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T); 
Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); 
PVE <- .6
rho <- .8
Pathways.Num <- 30
Pathways.Num.Selected <- 5
Pathways.SNPs <- 50
Pathways.SNPs.Interaction <- 200
Genome.Additive.Prop <- .5

#Additive SNPs
Data1.Additive <- Data1[,sample(1:ncol(Data1), ncol(Data1)*Genome.Additive.Prop)];
Data1.Additive.Betas <- rnorm(ncol(Data1.Additive),0,1);
Y.Additive <- Data1.Additive %*% Data1.Additive.Betas;
Data1.Additive.Betas <- Data1.Additive.Betas * c(sqrt((PVE * rho) / var(Y.Additive))); #Normalize in respect to PVE
Y.Additive <- Data1.Additive %*% Data1.Additive.Betas;

#Pathways & Genome Partners Make
Pathways <- c(); 
Pathways.Full <- sample(1:ncol(Data1), Pathways.Num * Pathways.SNPs); 
Count1 <- 1; for (i in 1:Pathways.Num) { 
	Pathways <- rbind(Pathways, Pathways.Full[Count1:(Count1+Pathways.SNPs-1)]); 
	Count1 <- Count1 + Pathways.SNPs - 1; 
}; 
##SNPs.Pathways <- c(); 
Genome.AntiPathway.SNPs <- c(); 
Genome.Total.SNPs <- 1:ncol(Data1); 
for (j in 1:Pathways.Num.Selected) { 
	Genome.AntiPathway.SNPs <- rbind(Genome.AntiPathway.SNPs, sample(Genome.Total.SNPs[! Genome.Total.SNPs %in% Pathways[j,]], Pathways.SNPs.Interaction)); 
};
#for (i in 1:5) { print(table(Pathways[i,] %in% Genome.AntiPathway.SNPs[i,])); }; #Data check

#Pathway Epistasis
Y.Epistasis <- 0;
Data1.Epistasis <- list();
Data1.Epistasis.Alphas <- c();
for (k in 1:Pathways.Num.Selected) {
	Data1.Epistasis.temp1 <- c();
	Data1.Epistasis.Alphas.temp1 <- c();
	for (l in 1:Pathways.SNPs) {
		for (m in 1:Pathways.SNPs.Interaction) {
			Epistasis1 <- (Data1[,Pathways[k,l]] * Data1[,Genome.AntiPathway.SNPs[k,m]]);
			Alpha1 <- rnorm(1,0,1);
			Y.Epistasis <- Y.Epistasis + Alpha1 * Epistasis1; 
			Data1.Epistasis.temp1 <- cbind(Data1.Epistasis.temp1, Epistasis1);
			Data1.Epistasis.Alphas.temp1 <- c(Data1.Epistasis.Alphas.temp1, Alpha1);
		};
	};
	Data1.Epistasis[[k]] <- Data1.Epistasis.temp1;
	Data1.Epistasis.Alphas <- cbind(Data1.Epistasis.Alphas, Data1.Epistasis.Alphas.temp1);
}; 
Epistasis.PVE.Rescale <- sqrt((PVE * (1-rho)) / var(Y.Epistasis));
Y.Epistasis <- 0;
Count2 <- 1; for (k in 1:Pathways.Num.Selected) { for (l in 1:Pathways.SNPs) { for (m in 1:Pathways.SNPs.Interaction) {
			Y.Epistasis <- Y.Epistasis + Data1.Epistasis.Alphas[Count2] * Data1.Epistasis[Count2] * Epistasis.PVE.Rescale; #Normalize in respect to PVE 
			Count2 <- Count2 + 1;
}; }; }; 

#Error
PVE.Error <- (1 - PVE) * (var(Y.Additive + Y.Epistasis) / PVE);
Y.Error <- rnorm(nrow(Data1),0,sqrt(PVE.Error));

#Final Phenotype
Y.Final <- Y.Additive + Y.Epistasis + Y.Error

#Pathway Formatting
Pathways.Edits <- apply(Pathways, 1, function(x) { return(paste(x, collapse=",")); }); Pathways.Edits <- cbind(1:length(Pathways.Edits), Pathways.Edits); Pathways.Edits <- cbind(rep("Pathway", nrow(Pathways.Edits)), Pathways.Edits); Pathways.Edits.2 <- apply(Pathways.Edits[,1:2], 1, function(x) { return(paste(x, collapse="")); }); Pathways.Edits <- cbind(Pathways.Edits.2, Pathways.Edits[,3]); 

#PVE Checks
PVE.Check.Linear <- var(Y.Additive) / var(Y.Final)
PVE.Check.Epistasis <- var(Y.Epistasis) / var(Y.Final)
PVE.Check.Error <- var(Y.Error) / var(Y.Final)
PVE.Check.Pathways <- c(); for (i in 1:length(Data1.Epistasis)) { PVE.Check.Pathways <- c(PVE.Check.Pathways, var(Data1.Epistasis[[i]] %*% Data1.Epistasis.Alphas[,i]) / var(Y.Final)); };
print(c(PVE.Check.Linear, PVE.Check.Epistasis, PVE.Check.Error, PVE.Check.Pathways));

#Output writing
write.table(Y.Final, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Pathways, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(Pathways.Edits, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.Edits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(Genome.AntiPathway.SNPs, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Genome_AntiPathway_SNPs.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#"

##Y <- rnorm(nrow(Data1)); 
##write.table(Y, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.Orig.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
##write.table(SNPs.Pathways, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(Y, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.Orig.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
#write.table(Y.new, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
#write.table(Pathways, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(Pathways.Edits, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(SNPs.Pathways, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(SNPs.Genome, file="/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#correction_factor = np.sqrt(self.pve*(1.0-self.rho)/np.var(self.y_pathway))
#        eps = (1.0-self.pve)*np.var(self.y_additive+self.y_pathway)/self.pve
#        self.y_err = np.random.normal(0,np.sqrt(eps),size=self.X.shape[0])
#        self.y = self.y_additive + self.y_pathway + self.y_err + self.y_pcs
#for (n in (Pathways.Num.Selected+1):Pathways.Num) {
#};
#	self.linear_pve = np.var(np.dot(self.X_additive,self.beta))/np.var(self.y)
#        self.epistatic_pve = np.var(self.y_pathway)/np.var(self.y)
#        self.pcs_pve = np.var(self.y_pcs)/np.var(self.y)
#        print("additive pve: ",self.linear_pve)
#        print("pathway pve: ",self.epistatic_pve)
#        print("pc pve: ",self.pcs_pve)
#        print("PVE per pathway")
#        for idx,alpha_i in enumerate(self.alpha):
#            print(np.var(np.dot(self.get_W(idx),alpha_i))/np.var(self.y))
#        print("error pve: ",np.var(self.y_err)/np.var(self.y))

##R -q -e"
#library("devtools"); devtools::load_all();
#set.seed(582724); 
#X <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T);
#Y <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", header=F);
#Pathways <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", header=F);
#SNPs.Pathways <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Pathways.txt", header=F);
#SNPs.Genome <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/SimData.SNPs_Genome.txt", header=F);
#Y2 <- cbind(Y,Y);
#Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
#library("Rcpp"); library("RcppArmadillo"); library("RcppParallel"); library("doParallel");
#sourceCpp("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/src/MAPITR.cpp");

#X.mean <- apply(X, 2, mean); X.sd <- apply(X, 2, sd); X.adj <- t((t(X)-X.mean)/X.sd); 
#summary(lm(as.matrix(Y) ~ X.adj[,1252] + X.adj[,956] + X.adj[,1252] * X.adj[,956]))

#MAPITRBase(as.matrix(Y2), as.matrix(t(X)), Pathways.Full[1:2], cores=1)
#MAPITRBaseTest(as.matrix(Y), as.matrix(t(X)), Pathways.Full[1], cores=1)
#MAPITRBaseTest(as.matrix(Y), Pathways.Full[1], cores=1)
#Test2 <- 

library("devtools"); devtools::load_all();
library("Rcpp"); library("RcppArmadillo"); library("RcppParallel"); library("doParallel"); library("CompQuadForm");
sourceCpp("/home/mturchin20/TempStuff/MAPITR/src/MAPITR.cpp")
set.seed(582724); 
X <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T);
Y <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.txt", header=F);
Pathways <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.Edits.txt", header=F);
SNPs.Pathways <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Pathways.txt", header=F);
SNPs.Genome <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Genome.txt", header=F);
Y2 <- cbind(Y,Y);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 

##MAPITRmain <- function (Phenotype, Genotypes, Pathways, GRM_Grand = NULL, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) 
#MAPITR.Results <- MAPITRmain(Y, X, Pathways);

Y30 <- matrix(unlist(rep(Y, 30)), ncol=30)
Results.temp2 <- MAPITRBase(Y30, t(X), Pathways.Full, cores=1)
#MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp2$Est, MAPITRoutput.temp2$Eigenvalues)
Results.temp2.pValues <- GetMAPITRpValues(Results.temp2$Est, Results.temp2$Eigenvalues)
Results.temp2.pValues


#30 pathways
#4 true
#each have 200 SNPs 
#10 real SNPs each interacting with 10 other SNPs in the genome

~~~
#20200627
> Data1.Additive.Betas <- rnorm(ncol(Data1.Additive),0,1);
> length(Data1.Additive.Betas)
[1] 5000
> Y.Additive <- Data1.Additive %*% Data1.Additive.Betas;
> length(Y.Additive)
[1] 2000
> sqrt((PVE * rho) / var(Y.Additive))
           [,1]
[1,] 0.00950565
> var(Y.Additive)
         [,1]
[1,] 5312.239
> Data1.Additive.Betas <- Data1.Additive.Betas * sqrt((PVE * rho) / var(Y.Additive)); #Normalize in respect to PVE
Warning message:
In Data1.Additive.Betas * sqrt((PVE * rho)/var(Y.Additive)) :
  Recycling array of length 1 in vector-array arithmetic is deprecated.
  Use c() or as.vector() instead.
  
> Y.Additive <- Data1.Additive %*% Data1.Additive.Betas;
> head(Data1.Additive.Betas * 2)
[1]  0.009188504 -0.002522482  0.036276107  0.001369265  0.023999342
[6] -0.031937548
> head(Data1.Additive.Betas * sqrt((PVE * rho) / var(Y.Additive)))
[1]  0.0045942520 -0.0012612410  0.0181380535  0.0006846324  0.0119996711
[6] -0.0159687741
Warning message:
In Data1.Additive.Betas * sqrt((PVE * rho)/var(Y.Additive)) :
  Recycling array of length 1 in vector-array arithmetic is deprecated.
  Use c() or as.vector() instead.

> sqrt((PVE * rho)/var(Y.Additive))
     [,1]
[1,]    1
> c(sqrt((PVE * rho)/var(Y.Additive)))
[1] 1
> var(Y.Additive)
     [,1]
[1,] 0.48
> (PVE * rho)
[1] 0.48
> dim(Pathways)
[1] 30 50
> Pathways.Num * Pathways.SNPs
[1] 1500
> length(Pathways.Full
+ )
[1] 1500
> Pathways[1:2,]
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
[1,] 9403 2298 1900 3218 4905 2302 2496  159 8667  1013  8392  1455  8661  6143
[2,] 3147 5575 8072  771 5471 2837 3894 3440 4016  6048  7685  6329  3831  9532
     [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
[1,]  5277  4125  7428  5330  7290  5308  8919  2469  3400   960  4391  6969
[2,]  5204  8701  2940  9810  5943  4329  6555  1665   988  7866  9970   395
     [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
[1,]  1465   527  5872  9641  3455  2762  8772  2986  3586  5430  1462  6739
[2,]  1690  7547  5339  5352  3382  3125  5370  4953  2828  4128  9407  1967
     [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50]
[1,]  6510  6938  1482  5130  6179  3499   327   193  1427  6843  7523  3147
[2,]  6262  1755   318  5552  9282  3870  9750  3779   270  7090  7608  4525
> Pathways.Full[1:51]
 [1] 9403 2298 1900 3218 4905 2302 2496  159 8667 1013 8392 1455 8661 6143 5277
[16] 4125 7428 5330 7290 5308 8919 2469 3400  960 4391 6969 1465  527 5872 9641
[31] 3455 2762 8772 2986 3586 5430 1462 6739 6510 6938 1482 5130 6179 3499  327
[46]  193 1427 6843 7523 3147 5575
> Pathways.Full[1:55]
 [1] 9403 2298 1900 3218 4905 2302 2496  159 8667 1013 8392 1455 8661 6143 5277
[16] 4125 7428 5330 7290 5308 8919 2469 3400  960 4391 6969 1465  527 5872 9641
[31] 3455 2762 8772 2986 3586 5430 1462 6739 6510 6938 1482 5130 6179 3499  327
[46]  193 1427 6843 7523 3147 5575 8072  771 5471 2837
> for (i in 1:5) { print(table(Pathways[i,] %in% Genome.AntiPathway.SNPs[i,])); };

FALSE
   50

FALSE
   50

FALSE
   50

FALSE
   50

FALSE
   50


~~~

