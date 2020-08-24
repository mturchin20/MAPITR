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
#zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "set.seed(248573); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 20000); Rows <- sample(1:nrow(Data1), 4000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs2.gz
zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 5000 | perl -lane 'print join("\t", @F[0..99999]);' | R -q -e "set.seed(248573); Data1 <- read.table(file('stdin'), header=T); Columns <- sample(1:ncol(Data1), 20000); Rows <- sample(1:nrow(Data1), 4000); Data2 <- Data1[Rows,Columns]; write.table(Data2, file=\"\", quote=FALSE, row.names=FALSE, col.names=TRUE);" | grep -v ^\> | gzip > /users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs3.gz

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
##set.seed(973459); Data1 <- read.table("/Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T); 
##R -q -e "
##set.seed(973459); Data1 <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz", header=T); 
##Data1 <- Data1[sample(1:nrow(Data1), 2000), sample(1:ncol(Data1), 10000)];
##set.seed(973459); Data1 <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T); 
#set.seed(973459); Data1 <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.gz", header=T);
set.seed(973459); Data1 <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs3.gz", header=T);
Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); 
PVE <- .8
rho <- .2
Pathways.Num <- 2
#Pathways.Num.Selected <- 2
Pathways.Num.Selected <- 1
Pathways.SNPs <- 100
Pathways.SNPs.Interaction <- 1000
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
Data1.Epistasis <- list(); Data1.Epistasis.Pathway.SNPs.Check <- list(); Data1.Epistasis.Genome.SNPs.Check <- list();
Data1.Epistasis.Alphas <- c();
for (k in 1:Pathways.Num.Selected) {
	Data1.Epistasis.temp1 <- c(); Data1.Epistasis.Pathway.SNPs.Check.temp1 <- c(); Data1.Epistasis.Genome.SNPs.Check.temp1 <- c();
	for (l in 1:Pathways.SNPs) {
		print(c(k,l));
		Data1.Epistasis.Pathway.SNPs <- c();
		Data1.Epistasis.Genome.SNPs <- c(); 
		for (m in 1:Pathways.SNPs.Interaction) {
			Data1.Epistasis.Pathway.SNPs <- cbind(Data1.Epistasis.Pathway.SNPs, Data1[,Pathways[k,l]]);
			Data1.Epistasis.Genome.SNPs <- cbind(Data1.Epistasis.Genome.SNPs, Data1[,Genome.AntiPathway.SNPs[k,m]]);
		};
		Epistasis1 <- Data1.Epistasis.Pathway.SNPs * Data1.Epistasis.Genome.SNPs;
		Data1.Epistasis.temp1 <- cbind(Data1.Epistasis.temp1, Epistasis1); Data1.Epistasis.Pathway.SNPs.Check.temp1 <- cbind(Data1.Epistasis.Pathway.SNPs.Check.temp1, Data1.Epistasis.Pathway.SNPs); Data1.Epistasis.Genome.SNPs.Check.temp1 <- cbind(Data1.Epistasis.Genome.SNPs.Check.temp1, Data1.Epistasis.Genome.SNPs); 
	};
	Data1.Epistasis.Alphas.temp1 <- rnorm(Pathways.SNPs*Pathways.SNPs.Interaction,0,1);
	Y.Epistasis <- Y.Epistasis + (Data1.Epistasis.temp1 %*% Data1.Epistasis.Alphas.temp1);
	Data1.Epistasis[[k]] <- Data1.Epistasis.temp1; Data1.Epistasis.Pathway.SNPs.Check[[k]] <- Data1.Epistasis.Pathway.SNPs.Check.temp1; Data1.Epistasis.Genome.SNPs.Check[[k]] <- Data1.Epistasis.Genome.SNPs.Check.temp1;
	Data1.Epistasis.Alphas <- cbind(Data1.Epistasis.Alphas, Data1.Epistasis.Alphas.temp1);
}; 
Epistasis.PVE.Rescale <- sqrt((PVE * (1-rho)) / var(Y.Epistasis));
#for (i in 1:length(Data1.Epistasis)) { print(dim(Data1.Epistasis[[i]])); }; print(dim(Data1.Epistasis.Alphas)); #Data check
#Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,198:207]; Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,398:407]; Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,598:607]; Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,798:807]; #Data check
#Data1.Epistasis.Pathway.SNPs.Check[[3]][1:5,198:207]; Data1.Epistasis.Pathway.SNPs.Check[[3]][1:5,398:407]; Data1.Epistasis.Pathway.SNPs.Check[[3]][1:5,598:607]; Data1.Epistasis.Pathway.SNPs.Check[[3]][1:5,798:807]; #Data check
#Data1.Epistasis.Genome.SNPs.Check[[1]][1:5,c(1:5,200:205)]; Data1.Epistasis.Genome.SNPs.Check[[1]][1:5,c(6:10,406:410)]; Data1.Epistasis.Genome.SNPs.Check[[1]][1:5,c(11:15,611:615)]; Data1.Epistasis.Genome.SNPs.Check[[1]][1:5,c(16:20,816:820)]; #Data check
#Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(1:5,200:205)]; Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(6:10,406:410)]; Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(11:15,611:615)]; Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(16:20,816:820)]; #Data check
Y.Epistasis <- 0;
for (k in 1:Pathways.Num.Selected) { 
	Y.Epistasis <- Y.Epistasis + ((Data1.Epistasis[[k]] %*% Data1.Epistasis.Alphas[,k]) * c(Epistasis.PVE.Rescale));
}; 

#Error
PVE.Error <- (1 - PVE) * (var(Y.Additive + Y.Epistasis) / PVE);
Y.Error <- rnorm(nrow(Data1),0,sqrt(PVE.Error));

#Final Phenotype
Y.Final <- Y.Additive + Y.Epistasis + Y.Error

#Pathway Formatting
Pathways.Edits <- apply(Pathways, 1, function(x) { return(paste(x, collapse=",")); }); Pathways.Edits <- cbind(1:length(Pathways.Edits), Pathways.Edits); Pathways.Edits <- cbind(rep("Pathway", nrow(Pathways.Edits)), Pathways.Edits); Pathways.Edits.2 <- apply(Pathways.Edits[,1:2], 1, function(x) { return(paste(x, collapse="")); }); Pathways.Edits <- cbind(Pathways.Edits.2, Pathways.Edits[,3]); 

#PVE Checks
PVE.Check.Additive <- var(Y.Additive) / var(Y.Final)
PVE.Check.Epistasis <- var(Y.Epistasis) / var(Y.Final)
PVE.Check.Error <- var(Y.Error) / var(Y.Final)
PVE.Check.Pathways <- c(); for (i in 1:length(Data1.Epistasis)) { PVE.Check.Pathways <- c(PVE.Check.Pathways, var(Data1.Epistasis[[i]] %*% Data1.Epistasis.Alphas[,i] * c(Epistasis.PVE.Rescale)) / var(Y.Final)); };
print(c(PVE, rho, PVE * rho, PVE - PVE * rho)); print(c(PVE.Check.Additive, PVE.Check.Epistasis, PVE.Check.Error)); print(PVE.Check.Pathways);

#Output writing
write.table(Y.Final, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Pathways, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(Pathways.Edits, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
write.table(Genome.AntiPathway.SNPs, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Genome_AntiPathway_SNPs.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
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
#write.table(Y.Final, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
#write.table(Pathways, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(Pathways.Edits, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Pathways.Edits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
#write.table(Genome.AntiPathway.SNPs, file="/home/mturchin20/TempStuff/MAPITR/SimData/SimData.Genome_AntiPathway_SNPs.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
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
#			Alpha1 <- rnorm(1,0,1);
#			Data1.Epistasis.Alphas.temp1 <- c(Data1.Epistasis.Alphas.temp1, Alpha1);
#			Y.Epistasis <- Y.Epistasis + Alpha1 * Epistasis1; 
#	Data1.Epistasis.Alphas.temp1 <- c();
#			Epistasis1 <- (Data1[,Pathways[k,l]] * Data1[,Genome.AntiPathway.SNPs[k,m]]);
#			Data1.Epistasis.temp1 <- cbind(Data1.Epistasis.temp1, Epistasis1);
#	Epistasis1 <- Data1.Epistasis.Pathway.SNPs * Data1.Epistasis.Genome.SNPs;
#	Data1.Epistasis[[k]] <- Epistasis1;
#Count2 <- 1; for (k in 1:Pathways.Num.Selected) { for (l in 1:Pathways.SNPs) { for (m in 1:Pathways.SNPs.Interaction) {
#			Y.Epistasis <- Y.Epistasis + Data1.Epistasis.Alphas[Count2] * Data1.Epistasis[Count2] * Epistasis.PVE.Rescale; #Normalize in respect to PVE 
#			Count2 <- Count2 + 1;
#}; }; }; 







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
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR/src/MAPITR.cpp")
set.seed(582724); 
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs3.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pheno.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Pathways.Edits.txt", header=F);
Genome.AntiPathway.SNPs <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData.Genome_AntiPathway_SNPs.txt", header=F);
#Y2 <- cbind(Y,Y);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
X.mean <- apply(X, 2, mean); X.sd <- apply(X, 2, sd); X.adj <- t((t(X)-X.mean)/X.sd); 

##MAPITRmain <- function (Phenotype, Genotypes, Pathways, GRM_Grand = NULL, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) 
#MAPITR.Results <- MAPITRmain(Y, X, Pathways);

#Y30 <- matrix(unlist(rep(Y, 30)), ncol=30)
Y30 <- c();
for (i in 1:nrow(Pathways)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(Y) ~ as.matrix(X.adj[,Pathways.Full[[i]]]) - 1))); };

Results.temp2 <- MAPITRBase(Y30, t(X.adj), Pathways.Full, cores=1)
#Results.temp2 <- MAPITRBaseOrig(Y30, t(X.adj), Pathways.Full, cores=1)
#MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp2$Est, MAPITRoutput.temp2$Eigenvalues)
Results.temp2.pValues <- GetMAPITRpValues(Results.temp2$Est, Results.temp2$Eigenvalues)
Results.temp2.pValues


#30 pathways
#4 true
#each have 200 SNPs 
#10 real SNPs each interacting with 10 other SNPs in the genome

#SNPs.Pathways <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Pathways.txt", header=F);
#SNPs.Genome <- read.table("/home/mturchin20/TempStuff/MAPITR/SimData/SimData.SNPs_Genome.txt", header=F);







#20200813
#mkdir /users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1
#scp -p /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Simulations/Simulations.zip mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/.
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_snp_list.RData")
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_ids.RData")
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/chromosome16_snps.RData")
#write.table(as.matrix(gene_snp_list), file="/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_snp_list.txt", quote=FALSE, row.names=FALSE)
#write.table(gene_ids, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_ids.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(chromosome16_snps, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/chromosome16_snps.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

mkdir /users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Temp1
cd /users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Temp1

#cp -p ../../MAPITR_temp1/Simulations/Code/InterPath.cpp ../../MAPITR_temp1/Simulations/Code/InterPath.edits1.cpp
#cp -p /users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.cpp /users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp

#R -q -e "
set.seed(11151990); library(doParallel); library(Rcpp); library(RcppArmadillo); library(RcppParallel); library(CompQuadForm); library(Matrix); library(MASS); library(truncnorm)
load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_snp_list.RData")
load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_ids.RData")
load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/chromosome16_snps.RData")

sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.cpp")

gene_list = list()

for(i in 1:ncol(gene_snp_list)){
  x = unlist(gene_snp_list[[i]])
  gene_list[[i]] = x[!is.na(x)]
  names(gene_list)[i] = colnames(gene_snp_list)[i]
}

X = chromosome16_snps; 
colnames(X) = sapply(colnames(X),function(x) unlist(strsplit(x,split = "_"))[1])
unique.snps = apply(X,2,function(x) length(unique(x)))
X = X[,unique.snps>1]; dim(X)
#X.lines <- apply(X, 2, function(x) { return(paste(x, collapse=",")); });
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

ind = nrow(X); nsnp = ncol(X)

### Make Sure that Gene List and Genotype SNPs Match Up ###
for(i in 1:length(gene_list)){
  x = gene_list[[i]]
  gene_list[[i]] = x[x%in%colnames(X)]
}

### Remove Duplicate SNPs ###
repeated_snps = unlist(gene_list)[duplicated(unlist(gene_list))]
for(i in 1:length(gene_list)){
  x = gene_list[[i]]
  gene_list[[i]] = x[which(x%in%repeated_snps==FALSE)]
}

### Remove Duplicate SNPs (based on geno) ###
X.lines <- apply(X, 2, function(x) { return(paste(x, collapse=",")); });
duplicated_snps <- names(X.lines[duplicated(X.lines)]) 
for(i in 1:length(gene_list)){
  x = gene_list[[i]]
  gene_list[[i]] = x[which(x%in%duplicated_snps==FALSE)]
}
#dim(X)
#X <- X[,!colnames(X)%in%names(X.lines[duplicated(X.lines)])]
#dim(X)
rm(X.lines)

### Remove Any Gene with One or Fewer SNPs ###
gene_list = gene_list[unlist(lapply(gene_list,function(x){length(x)}))>=2]

npthwy = length(gene_list) #Final Number of Genes/Pathways

### Define the Simulation Parameters ###
n.datasets = 1 #Total Number of Simulations
pve = 0.8; #Heritability of the trait
rho = 0.2; #Proportion of the heritability caused by additive effects {0.8, 0.5}

### Set Up Causal Pathways in Three Groups
ncausal1 = 5; ncausal2 = 50 #Pathways in each group
ncausal3 = 150-(ncausal1+ncausal2) #Total Pathways

### Create a list to save the final Results ###
pval_mat = matrix(nrow = npthwy,ncol = n.datasets); rownames(pval_mat) = names(gene_list)
G1_snps = matrix(nrow = ncausal1,ncol = n.datasets)
G2_snps = matrix(nrow = ncausal2,ncol = n.datasets)

### Run the Analysis ###
#for(j in 1:n.datasets){

  #Select Causal Pathways
  pthwy.ids = 1:npthwy
  s1=sample(pthwy.ids, ncausal1, replace=F)
  s2=sample(pthwy.ids[-s1], ncausal2, replace=F)
  s3=sample(pthwy.ids[-c(s1,s2)], ncausal3, replace=F)
  s1=seq(1,ncausal1,by=1)
  s2=seq(ncausal1+1,ncausal1+1+ncausal2,by=1)

  ### Simulate the Additive Effects ###
  snps = unlist(gene_list[c(s1,s2,s3)])
  Xmarginal = X[,snps]
  beta=rnorm(dim(Xmarginal)[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta

  ### Simulate Pairwise Interaction matrix ###
  Xepi = c(); b = c()
  causal_pthwys = gene_list[c(s1,s2,s3)]
  for(i in 1:ncausal1){
    snps = unlist(gene_list[s1[i]])
    for(k in 1:length(snps)){
      Xepi = cbind(Xepi,X[,snps[k]]*X[,unlist(gene_list[s2])])
    }
  }

  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta

  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve)/var(y_err))

  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_err

  ### Check dimensions ###
  dim(X); dim(y)

  ######################################################################################
  ######################################################################################
  ######################################################################################

  ### Rewrite the Pathway List to Give Column IDs Instead of SNP Names ###
  regions = list()
  for(i in 1:length(gene_list)){
    regions[[i]] = which(colnames(X)%in%gene_list[[i]])
    names(regions)[i] = names(gene_list)[i]
  }

  ######################################################################################
  ######################################################################################
  ######################################################################################

  ### Set the number of cores ###
  cores = detectCores()

#	save.image("20200813_temp1.RData")
#	save.image("20200813_temp2.RData")
#	save.image("20200813_temp3.RData")
#	save.image("20200813_temp4.RData") #diff pve/rho vals

set.seed(11151990); library(doParallel); library(Rcpp); library(RcppArmadillo); library(RcppParallel); library(CompQuadForm); library(Matrix); library(MASS); library(truncnorm)
	load("20200813_temp4.RData")
#	sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits1.cpp")
	sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp")
#	sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.cpp")
	regions <- regions[1:50]

#      b.cols(1,j.n_elem) = trans(X.rows(j-1));
#	X.t <- t(X);

Y30 <- c();
for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(y) ~ as.matrix(X[,regions[[i]]]) - 1))); };

  ### Run InterPath ###
  ptm <- proc.time() #Start clock
#  vc.mod = InterPath(t(X),y,regions,cores = cores)
  vc.mod = InterPath(t(X),Y30,regions,cores = cores)
  proc.time() - ptm #Stop clock

  ### Apply Davies Exact Method ###
  vc.ts = vc.mod$Est
  names(vc.ts) = names(regions)

  pvals = c()
  for(i in 1:length(vc.ts)){
    lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
    Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
    pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
    names(pvals)[i] = names(vc.ts[i])
  }
  pvals
 
  ### Find power for the first group of SNPs ###
  Pthwys_1 = names(regions)[s1]

  ### Find power for the second group of SNPs ###
  Pthwys_2 = names(regions)[s2]

  ######################################################################################
  ######################################################################################
  ######################################################################################

  ### Save Results ###
  pval_mat[,j] = pvals
  G1_snps[,j] = Pthwys_1
  G2_snps[,j] = Pthwys_2

  ### Report Status ###
  cat("Completed Dataset", j, "\n", sep = " ")
#}

#Save final Results
Final = list(pval_mat,G1_snps,G2_snps)

```
>   proc.time() - ptm #Stop clock
    user   system  elapsed 
1662.400   41.025 1703.461 
> pvals
        AARS         ABAT        ABCA3        ABCC1       ABCC11       ABCC12 
4.951986e-03 4.533479e-07 7.384488e-01 3.207512e-08 5.443513e-03 1.210351e-01 
       ABCC6        ACSF3        ACSM1       ACSM2A       ACSM2B        ACSM3 
3.571622e-02 1.828907e-03 3.114093e-01 3.372007e-01 7.752575e-01 6.586091e-01 
       ACSM5        ADAD2     ADAMTS18        ADAT1        ADCY7        ADCY9 
1.086930e-01 7.224191e-01 5.661781e-03 4.238015e-01 8.849778e-01 6.350021e-01 
      ADGRG1       ADGRG3       ADGRG5        AKTIP        ALDOA         ALG1 
1.511630e-02 9.641417e-04 6.753311e-01 4.773866e-02 8.054788e-01 8.621856e-01 
        AMFR      ANKRD11        ANKS3       ANKS4B        AP1G1        APOBR 
3.613904e-01 3.422696e-03 1.783861e-02 8.824898e-01 4.812328e-01 2.280406e-02 
      APOOP5         APRT         AQP8     ARHGAP17      ARL6IP1        ARMC5 
7.915088e-02 7.618808e-03 2.109255e-02 1.416500e-01 3.896289e-01 5.333194e-01 
     ATF7IP2        ATMIN       ATP2A1       ATP2C2     ATP6V0D1       ATXN1L 
8.343987e-01 9.534180e-01 8.267568e-01 4.671650e-05 5.826754e-01 9.069548e-01 
      ATXN2L        AXIN1       BAIAP3         BANP         BBS2        BCAR1 
5.068331e-01 6.661920e-01 1.898870e-01 2.785836e-02 3.677133e-01 4.780869e-01 
       BCL7C         BCO1 
3.490812e-02 4.476724e-01 
...
>         sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp")
>         regions <- regions[1:50]
>
>   ### Run InterPath ###
>   ptm <- proc.time() #Start clock
>   vc.mod = InterPath(t(X),y,regions,cores = cores)
>   proc.time() - ptm #Stop clock
    user   system  elapsed 
1661.762   40.384 1702.187 
> 
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est  
>   names(vc.ts) = names(regions)
> 
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
        AARS         ABAT        ABCA3        ABCC1       ABCC11       ABCC12 
4.951986e-03 4.533479e-07 7.384488e-01 3.207512e-08 5.443513e-03 1.210351e-01 
       ABCC6        ACSF3        ACSM1       ACSM2A       ACSM2B        ACSM3 
3.571622e-02 1.828907e-03 3.114093e-01 3.372007e-01 7.752575e-01 6.586091e-01 
       ACSM5        ADAD2     ADAMTS18        ADAT1        ADCY7        ADCY9 
1.086930e-01 7.224191e-01 5.661781e-03 4.238015e-01 8.849778e-01 6.350021e-01 
      ADGRG1       ADGRG3       ADGRG5        AKTIP        ALDOA         ALG1 
1.511630e-02 9.641417e-04 6.753311e-01 4.773866e-02 8.054788e-01 8.621856e-01 
        AMFR      ANKRD11        ANKS3       ANKS4B        AP1G1        APOBR 
3.613904e-01 3.422696e-03 1.783861e-02 8.824898e-01 4.812328e-01 2.280406e-02 
      APOOP5         APRT         AQP8     ARHGAP17      ARL6IP1        ARMC5 
7.915088e-02 7.618808e-03 2.109255e-02 1.416500e-01 3.896289e-01 5.333194e-01 
     ATF7IP2        ATMIN       ATP2A1       ATP2C2     ATP6V0D1       ATXN1L 
8.343987e-01 9.534180e-01 8.267568e-01 4.671650e-05 5.826754e-01 9.069548e-01 
      ATXN2L        AXIN1       BAIAP3         BANP         BBS2        BCAR1 
5.068331e-01 6.661920e-01 1.898870e-01 2.785836e-02 3.677133e-01 4.780869e-01 
       BCL7C         BCO1 
3.490812e-02 4.476724e-01 
#no covars/corrct
> #  vc.mod = InterPath(t(X),matrix(Y30, ncol=length(regions), byrow=FALSE),regions,cores = cores)
>   proc.time() - ptm #Stop clock
    user   system  elapsed
1474.966   42.586 1517.611
>
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est
>   names(vc.ts) = names(regions)
>
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
        AARS         ABAT        ABCA3        ABCC1       ABCC11       ABCC12
9.075350e-03 1.332589e-07 7.269296e-01 8.315985e-09 8.904773e-04 7.395582e-02
       ABCC6        ACSF3        ACSM1       ACSM2A       ACSM2B        ACSM3
4.724206e-02 1.221150e-03 3.537905e-01 3.758361e-01 8.866383e-01 7.561057e-01
       ACSM5        ADAD2     ADAMTS18        ADAT1        ADCY7        ADCY9
1.835457e-01 8.593086e-01 4.683992e-03 2.367160e-01 7.598773e-01 4.633264e-01
      ADGRG1       ADGRG3       ADGRG5        AKTIP        ALDOA         ALG1
8.094992e-03 7.423545e-04 5.489022e-01 4.195427e-02 6.875581e-01 9.103729e-01
        AMFR      ANKRD11        ANKS3       ANKS4B        AP1G1        APOBR
3.975499e-01 9.669478e-04 1.572000e-02 9.477620e-01 4.750039e-01 2.166131e-02
      APOOP5         APRT         AQP8     ARHGAP17      ARL6IP1        ARMC5
3.922667e-02 2.086286e-02 2.057187e-02 1.623535e-01 3.846440e-01 2.715451e-01
     ATF7IP2        ATMIN       ATP2A1       ATP2C2     ATP6V0D1       ATXN1L
9.739844e-01 9.879812e-01 8.828891e-01 9.662260e-05 7.345508e-01 7.031548e-01
      ATXN2L        AXIN1       BAIAP3         BANP         BBS2        BCAR1
6.785574e-01 7.652838e-01 1.987019e-01 3.239428e-02 3.443410e-01 5.144689e-01
       BCL7C         BCO1
4.128038e-02 3.706294e-01
#with regress out before
>   pvals
        AARS         ABAT        ABCA3        ABCC1       ABCC11       ABCC12 
1.069194e-02 2.052235e-06 8.588981e-01 9.648794e-08 1.368308e-02 2.240554e-01 
       ABCC6        ACSF3        ACSM1       ACSM2A       ACSM2B        ACSM3 
5.389115e-02 6.202295e-03 2.422489e-01 3.721679e-01 6.374053e-01 7.814415e-01 
       ACSM5        ADAD2     ADAMTS18        ADAT1        ADCY7        ADCY9 
6.661380e-02 6.291496e-01 1.482319e-02 5.528941e-01 9.797160e-01 8.027934e-01 
      ADGRG1       ADGRG3       ADGRG5        AKTIP        ALDOA         ALG1 
2.129916e-02 1.458324e-03 7.566949e-01 6.699238e-02 8.698268e-01 9.038601e-01 
        AMFR      ANKRD11        ANKS3       ANKS4B        AP1G1        APOBR 
4.557393e-01 1.116736e-02 3.273128e-02 9.573315e-01 6.185752e-01 1.745093e-02 
      APOOP5         APRT         AQP8     ARHGAP17      ARL6IP1        ARMC5 
1.200981e-01 2.040583e-02 2.766836e-02 1.940421e-01 4.164473e-01 6.224320e-01 
     ATF7IP2        ATMIN       ATP2A1       ATP2C2     ATP6V0D1       ATXN1L 
6.005474e-01 8.230449e-01 8.874864e-01 1.022908e-04 9.108130e-01 9.659642e-01 
      ATXN2L        AXIN1       BAIAP3         BANP         BBS2        BCAR1 
4.662691e-01 5.283394e-01 2.597461e-01 3.938402e-02 4.056004e-01 7.111617e-01 
       BCL7C         BCO1 
5.787721e-02 3.685350e-01 



```
 








set.seed(686234); library(doParallel); library(Rcpp); library(RcppArmadillo); library(RcppParallel); library(CompQuadForm); library(Matrix); library(MASS); library(truncnorm)
Data1 <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs3.gz", header=T);
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp")

#MAFs <- colSums(Data1)/(2*nrow(Data1));
#Data2 <- c(); for (i in 1:ncol(Data1)) { Data2 <- cbind(Data2, rbinom(nrow(Data1), 2, MAFs[i])); };
Data3 <- apply(Data1, 2, function(x) { MAF <- sum(x)/(2*length(x)); return(rbinom(length(x), 2, MAF)); }); 
#Data2 <- apply(Data1, 2, function(x) { return(sample(x)); });
X = Data3; 
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

ind = nrow(X); nsnp = ncol(X)

### Define the Simulation Parameters ###
n.datasets = 1 #Total Number of Simulations
pve = 0.6; #Heritability of the trait
rho = 0.6; #Proportion of the heritability caused by additive effects {0.8, 0.5}

### Set Up Causal SNPs
n.pathways = 2 #Number of Pathways
ncausal1a = 50; ncausal1b = 750 #
ncausal2a = 100; ncausal2b = 1000 #

  #Select Causal Pathways
  snp.ids = 1:ncol(X)
  s1a=sample(snp.ids, ncausal1a, replace=F)
  s1b=sample(snp.ids[-s1a], ncausal1b, replace=F)
#  s2a=sample(snp.ids, ncausal2a, replace=F)
#  s2b=sample(snp.ids[-s2a], ncausal2b, replace=F)

  regions <- list(); regions[[1]] <- s1a;
  for (i in 2:5) {
	regions[[i]] <- sample(snp.ids, 50, replace=F);
  }

  #Additive Effects
#  snps = c(s1a,s1b,s2a,s2b); 
  snps1 = unlist(regions[c(1,2,3,4,5)])
#  snps1 = c(s1a,s1b); 
  Xmarginal = X[,snps1]
  beta=rnorm(dim(Xmarginal)[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
 
  ### Simulate Pairwise Interaction matrix ###
  Xepi = c(); b = c()
  snps2 = s1a; 
  for(k in 1:length(snps2)){
    print(k)
    Xepi = cbind(Xepi,X[,snps2[k]]*X[,s1b])
  }

  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta

  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve)/var(y_err))

  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_err

  ### Check dimensions ###
  dim(X); dim(y)

#Pathway Formatting
Pathways <- matrix(unlist(regions), nrow=5, byrow=5)
Pathways.Edits <- apply(Pathways, 1, function(x) { return(paste(x, collapse=",")); }); Pathways.Edits <- cbind(1:length(Pathways.Edits), Pathways.Edits); Pathways.Edits <- cbind(rep("Pathway", nrow(Pathways.Edits)), Pathways.Edits); Pathways.Edits.2 <- apply(Pathways.Edits[,1:2], 1, function(x) { return(paste(x, collapse="")); }); Pathways.Edits <- cbind(Pathways.Edits.2, Pathways.Edits[,3]); 

cores = detectCores()
Y30 <- c();
for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(y) ~ as.matrix(X[,regions[[i]]]) - 1))); };

  ptm <- proc.time() #Start clock
  vc.mod = InterPath(t(X),Y30,regions,cores = cores)
  proc.time() - ptm #Stop clock
  
  ### Apply Davies Exact Method ###
  vc.ts = vc.mod$Est
  names(vc.ts) = names(regions)

  pvals = c()
  for(i in 1:length(vc.ts)){
    lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
    Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
    pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
    names(pvals)[i] = names(vc.ts[i])
  }
  pvals

#Output writing
write.table(Data3, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Genotypes.txt", quote=FALSE, row.names=FALSE, col.names=TRUE);
write.table(y, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Phenotype.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Pathways.Edits, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
system("gzip -f /users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Genotypes.txt");
#"

library(doParallel); library(Rcpp); library(RcppArmadillo); library(RcppParallel); library(CompQuadForm); library(Matrix); library(MASS); library(truncnorm)
#sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp")
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR/src/MAPITR.cpp")
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Genotypes.txt.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Phenotype.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Pathways.txt", header=F);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
regions <- Pathways.Full;
regions2 <- lapply(regions, function(x) { return(x-1);});
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
cores = detectCores()

Y30 <- c();
for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(Y) ~ as.matrix(X[,regions[[i]]]) - 1))); };

  ptm <- proc.time() #Start clock
#  vc.mod = InterPath(t(X),Y30,regions,cores = cores)
  vc.mod = MAPITRBase(t(X),Y30,regions2,cores = cores)
  proc.time() - ptm #Stop clock
  
  ### Apply Davies Exact Method ###
  vc.ts = vc.mod$Est
  names(vc.ts) = names(regions)

  pvals = c()
  for(i in 1:length(vc.ts)){
    lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
    Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
    pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
    names(pvals)[i] = names(vc.ts[i])
  }
  pvals

#Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); 
#
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_snp_list.RData")
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/gene_ids.RData")
#load("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Data/chromosome16_snps.RData")
#
#gene_list = list()
#
#for(i in 1:ncol(gene_snp_list)){
#  x = unlist(gene_snp_list[[i]])
#  gene_list[[i]] = x[!is.na(x)]
#  names(gene_list)[i] = colnames(gene_snp_list)[i]
#}
#  regions <- list(); regions[[1]] <- s1a;
#
#	save.image("20200814_vs2_temp1.RData") #diff pve/rho vals

#From: https://github.com/r-spatialecology/landscapemetrics/issues/52, https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define, https://stackoverflow.com/questions/21944695/rcpparmadillo-and-arma-namespace, https://github.com/RcppCore/RcppArmadillo/pull/88, https://github.com/gvegayon/arma64bit, https://cran.r-project.org/web/packages/RcppArmadillo/RcppArmadillo.pdf
#library(Rcpp); library(RcppArmadillo); library(RcppParallel); library(CompQuadForm); library(Matrix); library(MASS); library(truncnorm)
library("devtools"); devtools::load_all(); 
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Genotypes.txt.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Phenotype.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Pathways.txt", header=F);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
regions <- Pathways.Full;
regions2 <- lapply(regions, function(x) { return(x-1);});
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
cores = parallel::detectCores()

Y30 <- c();
for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(Y) ~ as.matrix(X[,regions[[i]]]) - 1))); };

  ptm <- proc.time() #Start clock
  vc.mod = MAPITRBase(t(X),Y30,regions,cores = cores)
  proc.time() - ptm #Stop clock
  
  ### Apply Davies Exact Method ###
  vc.ts = vc.mod$Est
  names(vc.ts) = names(regions)

  pvals = c()
  for(i in 1:length(vc.ts)){
    lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
    Davies_Method = CompQuadForm::davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
    pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
    names(pvals)[i] = names(vc.ts[i])
  }
  pvals

#Data1.mean <- apply(Data1, 2, mean); Data1.sd <- apply(Data1, 2, sd); Data1 <- t((t(Data1)-Data1.mean)/Data1.sd); 










#Going from main MAPITR functions

#From: https://stackoverflow.com/questions/733926/openmp-coding-warning-ignoring-pragma-omp-parallel
library("devtools"); devtools::load_all(); 
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Genotypes.txt.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Phenotype.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Pathways.txt", header=F);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
regions <- Pathways.Full;

#Output1 <- MAPITR(X, Y, regions); 
Output1 <- MAPITRmain(X, Y, Pathways); 

Genotypes <- X
Phenotype <- Y
#Pathways <- Pathways
GRM_Grand = NULL; GRM_Pathway = NULL; Covariates = NULL; CenterStandardize = TRUE; RegressPhenotypes = TRUE; PrintProgress = FALSE

        MAPITRprocessing <- list()
        MAPITRoutput <- list()
        MAPITRoutput$LogFile <- c()
        MAPITRoutput$pValues <- NULL
        MAPITRoutput$PVE <- NULL
        cores = parallel::detectCores()

        MAPITRoutput$LogFile <- rbind(MAPITRoutput$LogFile, paste(format(Sys.time()), " -- beginning MAPITR.", sep=""))
        if (PrintProgress == TRUE) {
                write(paste(format(Sys.time()), " -- beginning MAPITR.", sep=""), stderr())
        }

        #Data Checks
        MAPITRoutput$LogFile <- DataChecks(Phenotype, Genotypes, Pathways, Covariates, MAPITRoutput$LogFile)

        #Preprocessing Data
        MAPITRoutput.temp1 <- PreprocessData(Phenotype, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, MAPITRoutput$LogFile)
        PhenotypeMatrix <- MAPITRoutput.temp1$PhenotypeMatrix
        Genotypes <- MAPITRoutput.temp1$Genotypes
        Pathways.Full <- MAPITRoutput.temp1$Pathways.Full
        MAPITRoutput$LogFile <- MAPITRoutput.temp1$LogFile
        rm(MAPITRoutput.temp1)
        print("yaya2")

#	save.image("20200820_temp1.RData")
library("devtools"); 
devtools::load_all(); 
	load("20200820_temp1.RData")

	RunMAPITR.Base.Output.temp1 <- MAPITRBase(t(as.matrix(Genotypes)),as.matrix(PhenotypeMatrix),Pathways.Full,cores=cores)
#	RunMAPITR.Base.Output.temp2 <- MAPITRBase2(t(Genotypes),PhenotypeMatrix,Pathways.Full,cores=cores)
#	RunMAPITR.Base.Output.temp3 <- MAPITRBase3(t(Genotypes),PhenotypeMatrix,Pathways.Full,cores=cores)
#	RunMAPITR.Base.Output.temp4 <- MAPITRBase4(t(Genotypes))

  vc.ts = RunMAPITR.Base.Output.temp1$Est
  names(vc.ts) = names(Pathways.Full)
  pvals = c()
  for(i in 1:length(vc.ts)){
    lambda = sort(RunMAPITR.Base.Output.temp1$Eigenvalues[,i],decreasing = T)
    Davies_Method = CompQuadForm::davies(RunMAPITR.Base.Output.temp1$Est[i], lambda = lambda, acc=1e-8)
    pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
    names(pvals)[i] = names(vc.ts[i])
  }
  pvals





#Going from main MAPITR functions

#From: https://stackoverflow.com/questions/733926/openmp-coding-warning-ignoring-pragma-omp-parallel
library("devtools"); devtools::load_all(); 
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Genotypes.txt.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Phenotype.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Pathways.txt", header=F);
Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric); 
regions <- Pathways.Full;

#Output1 <- MAPITR(X, Y, regions); 
Output1 <- MAPITRmain(X, Y, Pathways); 
Output2 <- MAPITR(X, Y, Pathways); 

MAPITR_SimData_Genotypes <- X;
MAPITR_SimData_Phenotype <- Y;
MAPITR_SimData_Pathways <- Pathways;

#20200823 NOTE -- moving to use 'TestData' as 'SimData' since 'SimData' is too large, throws warnings in the devtools::check() process
#saveCompress1 <- "bzip2"
#save(MAPITR_SimData_Genotypes, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Genotypes.rda", compresss=saveCompress1);
#save(MAPITR_SimData_Phenotype, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Phenotype.rda");
#save(MAPITR_SimData_Pathways, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Pathways.rda");







# Trying to create a unit test 'sized' simulated dataset

set.seed(173463); 
library("devtools"); devtools::load_all(); 
Data1 <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/ukb_chrAll_v3.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Simulation.cutdwn.vs3.gz", header=T);

set.seed(173463); 
Data2 <- Data1[sample(1:nrow(Data1), 750), sample(1:ncol(Data1), 1000)];
#Data3 <- apply(Data2, 2, function(x) { MAF <- sum(x)/(2*length(x)); return(rbinom(length(x), 2, MAF)); }); 
#Data3 <- c(); for (i in 1:ncol(Data2)) { Data3 <- cbind(Data3, rbinom(nrow(Data2), 2, runif(1,.4,.6))); };
Data3 <- c(); for (i in 1:ncol(Data2)) { Data3 <- cbind(Data3, rbinom(nrow(Data2), 2, .5)); };
X = Data3; 
#X <- Data3[sample(1:nrow(Data3), 100), sample(1:ncol(Data3), 1000)];

Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

ind = nrow(X); nsnp = ncol(X)

### Define the Simulation Parameters ###
n.datasets = 1 #Total Number of Simulations
pve = 0.9; #Heritability of the trait
rho = 0.1; #Proportion of the heritability caused by additive effects {0.8, 0.5}

### Set Up Causal SNPs
n.pathways = 2 #Number of Pathways
ncausal1a = 75; ncausal1b = 250 #
ncausal2a = 100; ncausal2b = 1000 #

  #Select Causal Pathways
  snp.ids = 1:ncol(X)
  s1a=sample(snp.ids, ncausal1a, replace=F)
  s1b=sample(snp.ids[-s1a], ncausal1b, replace=F)

  regions <- list(); regions[[1]] <- s1a;
  for (i in 2:5) {
	regions[[i]] <- sample(snp.ids, 75, replace=F);
  }

  #Additive Effects
  snps1 = unlist(regions[c(1,2,3,4,5)])
  Xmarginal = X[,snps1]
  beta=rnorm(dim(Xmarginal)[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
 
  ### Simulate Pairwise Interaction matrix ###
  Xepi = c(); b = c()
  snps2 = s1a; 
  for(k in 1:length(snps2)){
    print(k)
    Xepi = cbind(Xepi,X[,snps2[k]]*X[,s1b])
  }

  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta

  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve)/var(y_err))

  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_err

  ### Check dimensions ###
  dim(X); dim(y)

#Pathway Formatting
Pathways <- matrix(unlist(regions), nrow=5, byrow=5)
Pathways.Edits <- apply(Pathways, 1, function(x) { return(paste(x, collapse=",")); }); Pathways.Edits <- cbind(1:length(Pathways.Edits), Pathways.Edits); Pathways.Edits <- cbind(rep("Pathway", nrow(Pathways.Edits)), Pathways.Edits); Pathways.Edits.2 <- apply(Pathways.Edits[,1:2], 1, function(x) { return(paste(x, collapse="")); }); Pathways.Edits <- cbind(Pathways.Edits.2, Pathways.Edits[,3]); 

Output1 <- MAPITRmain(X, y, Pathways.Edits); 
Output1$Results

#Output writing
write.table(Data3, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Genotypes.txt", quote=FALSE, row.names=FALSE, col.names=TRUE);
write.table(y, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Phenotype.txt", quote=FALSE, row.names=FALSE, col.names=FALSE); 
write.table(Pathways.Edits, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Pathways.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);  
system("gzip -f /users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Genotypes.txt");

library("devtools"); devtools::load_all(); 
X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Genotypes.txt.gz", header=T);
Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Phenotype.txt", header=F);
Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Pathways.txt", header=F);

#Output1 <- MAPITRmain(X, Y, Pathways); 
#Output1$Results
Output2 <- MAPITR(X, Y, Pathways); 
Output2$Results

MAPITR_TestData_Genotypes <- X;
MAPITR_TestData_Phenotype <- Y;
MAPITR_TestData_Pathways <- Pathways;
MAPITR_SimData_Genotypes <- MAPITR_TestData_Genotypes
MAPITR_SimData_Phenotype <- MAPITR_TestData_Phenotype
MAPITR_SimData_Pathways <- MAPITR_TestData_Pathways

saveCompress2 <- "xz"
save(MAPITR_TestData_Genotypes, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_TestData_Genotypes.rda", compress=saveCompress2);
save(MAPITR_TestData_Phenotype, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_TestData_Phenotype.rda");
save(MAPITR_TestData_Pathways, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_TestData_Pathways.rda");
save(MAPITR_SimData_Genotypes, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Genotypes.rda", compress=saveCompress2);
save(MAPITR_SimData_Phenotype, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Phenotype.rda");
save(MAPITR_SimData_Pathways, file="/users/mturchin/LabMisc/RamachandranLab/MAPITR/data/MAPITR_SimData_Pathways.rda");



# Sim data for unit tests

#From: https://stackoverflow.com/questions/54056594/cran-acceptable-way-of-linking-to-openmp-some-c-code-called-from-rcpp, https://github.com/r-lib/testthat/issues/361
set.seed(183474)
vals1 <- rnorm(10,0,1)
vals2 <- rbinom(10,2,.75)
vals3 <- vals1 + (runif(10,0,4) * vals2)
resids <- residuals(lm(vals3 ~ vals2))
summary(lm(vals3 ~ vals2))

library("devtools"); 
devtools::load_all(); 
devtools::test();
devtools::build();



# Final steps for building package, webpage, and misc work for eventual CRAN upload
#From: https://pkgdown.r-lib.org/, https://sahirbhatnagar.com/blog/2020/03/03/creating-a-website-for-your-r-package/, https://stackoverflow.com/questions/34585560/travis-ci-r-package-error-in-documentation
library("devtools"); 
library("pkgdown");
#Note -- R data object documention needs to be comlpete before the `install()` process can properly finish
devtools::install()
pkgdown::build_site()
devtools::build_vignettes()
devtools::document()
devtools::build()
devtools::check()

library("devtools"); 
devtools::check()

#add the following to '.travis.yml':
```
language: R
cache: packages
latex: false
r_packages:
 - doParallel
 - Rcpp
 - RcppArmadillo
 - RcppParallel
 - CompQuadForm
r_build_args: "--no-build-vignettes"
r_check_args: "--no-build-vignettes --as-cran" 
```

#couldn't get any of this to work, use travis ci build and do local build on mac
#cd /users/mturchin/Software
#wget https://sourceforge.net/projects/qpdf/files/qpdf/10.0.1/qpdf-10.0.1.tar.gz
#tar -xvzf qpdf-10.0.1.tar.gz
#cd qpdf-10.0.1
#./configure --prefix=/users/mturchin/local
#make -I/users/mturchin/local/include/
#make install
#
#wget https://sourceforge.net/projects/libjpeg-turbo/files/2.0.5/libjpeg-turbo-2.0.5.tar.gz
#tar -xvzf libjpeg-turbo-2.0.5.tar.gz
#cd libjpeg-turbo-2.0.5
##cmake /users/mturchin/Software/libjpeg-turbo-2.0.5
##cmake --prefix /users/mturchin/Software/libjpeg-turbo-2.0.5 --install
#cmake -DCMAKE_INSTALL_PREFIX=/users/mturchin/local /users/mturchin/Software/libjpeg-turbo-2.0.5
#make
#make install

W  checking compilation flags in Makevars ...
   Non-portable flags in variable 'PKG_CXXFLAGS':
     -Wall -fopenmp

 Found the following hidden files and directories:
    .travis.yml
  These were most likely included in error. See section 'Package
  structure' in the 'Writing R Extensions' manual.
  






```
#.7, .4, 100/1000
>   proc.time() - ptm #Stop clock
   user  system elapsed
176.399   4.990 181.396
>
>
>
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est
>   names(vc.ts) = names(regions)
>
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
[1] 7.466439e-09 1.523614e-01 4.545533e-01 5.005918e-01 6.321546e-01
#.6, .7, 50/750
>   proc.time() - ptm #Stop clock
   user  system elapsed
176.636   4.929 181.572
>
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est
>   names(vc.ts) = names(regions)
>
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
[1] 3.647303e-06 8.572264e-01 2.054071e-01 7.735247e-01 7.712807e-01
#.6, .6, 50/750 -- using input files from sim output
> sourceCpp("/users/mturchin/LabMisc/RamachandranLab/MAPITR_temp1/Simulations/Code/InterPath.edits2.cpp")
> X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Genotypes.txt.gz", header=T);
> Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Phenotype.txt", header=F);
> Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/SimData/SimData2.Pathways.txt", header=F);
> Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric);
> regions <- Pathways.Full;
> Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
> cores = detectCores()
>
> Y30 <- c();
> for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(Y) ~ as.matrix(X[,regions[[i]]]) - 1))); };
>
>   ptm <- proc.time() #Start clock
>   vc.mod = InterPath(t(X),Y30,regions,cores = cores)
>   proc.time() - ptm #Stop clock
   user  system elapsed
237.176   8.100 245.310
>
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est
>   names(vc.ts) = names(regions)
>
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
[1] 5.137639e-07 1.932025e-01 5.878403e-01 8.407904e-01 9.580596e-01
#From main MAPITR anaysis function directly
> library("devtools"); devtools::load_all();
Loading required package: usethis

Loading MAPITR

Attaching package: 'testthat'

The following object is masked from 'package:devtools':

    test_file

> X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Genotypes.txt.gz", header=T);
> Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Phenotype.txt", header=F);
> Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData2.Pathways.txt", header=F);
> Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric);
> regions <- Pathways.Full;
> regions2 <- lapply(regions, function(x) { return(x-1);});
> Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
> cores = parallel::detectCores()
>
> Y30 <- c();
> for (i in 1:length(regions)) { Y30 <- cbind(Y30, residuals(lm(as.matrix(Y) ~ as.matrix(X[,regions[[i]]]) - 1))); };
>
>   ptm <- proc.time() #Start clock
>   vc.mod = MAPITRBase(t(X),Y30,regions,cores = cores)
>   proc.time() - ptm #Stop clock
   user  system elapsed
187.646   3.237 190.915
>
>   ### Apply Davies Exact Method ###
>   vc.ts = vc.mod$Est
>   names(vc.ts) = names(regions)
>
>   pvals = c()
>   for(i in 1:length(vc.ts)){
+     lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
+     Davies_Method = CompQuadForm::davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
+     pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
+     names(pvals)[i] = names(vc.ts[i])
+   }
>   pvals
[1] 5.137640e-07 1.932025e-01 5.878403e-01 8.407904e-01 9.580596e-01
#unit testing version
> library("devtools"); devtools::load_all();
Loading required package: usethis
Loading MAPITR

Attaching package: 'testthat'

The following object is masked from 'package:devtools':

    test_file

> X <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Genotypes.txt.gz", header=T);
> Y <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Phenotype.txt", header=F);
> Pathways <- read.table("/users/mturchin/LabMisc/RamachandranLab/MAPITR/temp1/SimData/SimData3.Pathways.txt", header=F);
> 
> Output1 <- MAPITRmain(X, Y, Pathways);
> Output1$Results
     Pathways      pValues       Est       PVE
[1,]        1 0.0001933362 1.0761720 1.2417707
[2,]        2 0.4721117327 0.2013154 0.2303672
[3,]        3 0.4745897324 0.2039635 0.2295905
[4,]        4 0.4314853397 0.2232131 0.2545830
[5,]        5 0.5693362759 0.1596100 0.1858096
# Sim data for PreProcessing functionality unit tests
> set.seed(183474)
> vals1 <- rnorm(10,0,1)
> vals2 <- rbinom(10,2,.75)
> vals3 <- vals1 + (runif(10,0,4) * vals2)
> resids <- residuals(lm(vals3 ~ vals2))
> summary(lm(vals3 ~ vals2))

Call:
lm(formula = vals3 ~ vals2)

Residuals:
    Min      1Q  Median      3Q     Max
-3.0887 -0.9650  0.2395  0.8895  3.0852

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -0.8156     1.3002  -0.627   0.5480
vals2         2.8065     0.8573   3.274   0.0113 *
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 2.117 on 8 degrees of freedom
Multiple R-squared:  0.5726,    Adjusted R-squared:  0.5191
F-statistic: 10.72 on 1 and 8 DF,  p-value: 0.0113
# unit tests test
> library("devtools");
Loading required package: usethis
> devtools::load_all();
Loading MAPITR

Attaching package: 'testthat'

The following object is masked from 'package:devtools':

    test_file

> devtools::test();
Loading MAPITR
Testing MAPITR
v |  OK F W S | Context
v |  26       | Basic R Functionality
v |   5       | Tests for DataPreprocessing.R
v |   5       | Tests for MAPITR.R [2.0 s]

== Results =====================================================================
Duration: 2.1 s

OK:       36
Failed:   0
Warnings: 0
Skipped:  0

You are a coding rockstar!

```












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
> vals1 <- matrix(c(1,2,3,4), ncol=2)
> vals2 <- matrix(c(6,5,7,8), ncol=2)
> vals1
     [,1] [,2]
[1,]    1    3
[2,]    2    4
> vals2
     [,1] [,2]
[1,]    6    7
[2,]    5    8
> vals1 * vals2
     [,1] [,2]
[1,]    6   21
[2,]   10   32
> Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,198:207]; Data1.Epistasis.Pathway.SNPs.Check[[1]][1:5,398:407];  
           [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
[1,] -0.6466086 -0.6466086 -0.6466086 -0.3680773 -0.3680773 -0.3680773
[2,]  1.1722622  1.1722622  1.1722622  2.5531709  2.5531709  2.5531709
[3,] -0.6466086 -0.6466086 -0.6466086  2.5531709  2.5531709  2.5531709
[4,] -0.6466086 -0.6466086 -0.6466086 -0.3680773 -0.3680773 -0.3680773
[5,] -0.6466086 -0.6466086 -0.6466086 -0.3680773 -0.3680773 -0.3680773
           [,7]       [,8]       [,9]      [,10]
[1,] -0.3680773 -0.3680773 -0.3680773 -0.3680773
[2,]  2.5531709  2.5531709  2.5531709  2.5531709
[3,]  2.5531709  2.5531709  2.5531709  2.5531709
[4,] -0.3680773 -0.3680773 -0.3680773 -0.3680773
[5,] -0.3680773 -0.3680773 -0.3680773 -0.3680773
           [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
[1,] -0.3680773 -0.3680773 -0.3680773  1.2272088  1.2272088  1.2272088
[2,]  2.5531709  2.5531709  2.5531709 -0.6680944 -0.6680944 -0.6680944
[3,]  2.5531709  2.5531709  2.5531709 -0.6680944 -0.6680944 -0.6680944
[4,] -0.3680773 -0.3680773 -0.3680773 -0.6680944 -0.6680944 -0.6680944
[5,] -0.3680773 -0.3680773 -0.3680773 -0.6680944 -0.6680944 -0.6680944
           [,7]       [,8]       [,9]      [,10]
[1,]  1.2272088  1.2272088  1.2272088  1.2272088
[2,] -0.6680944 -0.6680944 -0.6680944 -0.6680944
[3,] -0.6680944 -0.6680944 -0.6680944 -0.6680944
[4,] -0.6680944 -0.6680944 -0.6680944 -0.6680944
[5,] -0.6680944 -0.6680944 -0.6680944 -0.6680944
> Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(1:5,200:205)]; Data1.Epistasis.Genome.SNPs.Check[[4]][1:5,c(6:10,406:410)]; 
          [,1]      [,2]     [,3]       [,4]       [,5]      [,6]      [,7]
[1,] -0.508133 1.6492337 -0.27114 -1.0218675 -1.2463003  1.393562 -0.508133
[2,] -0.508133 0.2194689 -0.27114  1.9836252 -1.2463003 -0.590156 -0.508133
[3,]  1.730338 0.2194689 -0.27114  0.4808788  0.1579817 -0.590156  1.730338
[4,] -0.508133 0.2194689 -0.27114  0.4808788  1.5622638  1.393562 -0.508133
[5,] -0.508133 0.2194689 -0.27114 -1.0218675  1.5622638 -0.590156 -0.508133
          [,8]     [,9]      [,10]      [,11]
[1,] 1.6492337 -0.27114 -1.0218675 -1.2463003
[2,] 0.2194689 -0.27114  1.9836252 -1.2463003
[3,] 0.2194689 -0.27114  0.4808788  0.1579817
[4,] 0.2194689 -0.27114  0.4808788  1.5622638
[5,] 0.2194689 -0.27114 -1.0218675  1.5622638
           [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
[1,] -0.5145001 -0.1903977 -0.1979011  0.244801 -0.2859267 -0.5145001
[2,] -0.5145001 -0.1903977 -0.1979011 -1.178461 -0.2859267 -0.5145001
[3,] -0.5145001 -0.1903977 -0.1979011  0.244801 -0.2859267 -0.5145001
[4,] -0.5145001 -0.1903977 -0.1979011  0.244801 -0.2859267 -0.5145001
[5,]  1.5897949 -0.1903977 -0.1979011  0.244801 -0.2859267  1.5897949
           [,7]       [,8]      [,9]      [,10]
[1,] -0.1903977 -0.1979011  0.244801 -0.2859267
[2,] -0.1903977 -0.1979011 -1.178461 -0.2859267
[3,] -0.1903977 -0.1979011  0.244801 -0.2859267
[4,] -0.1903977 -0.1979011  0.244801 -0.2859267
[5,] -0.1903977 -0.1979011  0.244801 -0.2859267
#20200813
> length(X.lines[duplicated(X.lines)])
[1] 104
> names(X.lines[duplicated(X.lines)])
  [1] "16:180529"   "16:188768"   "16:633354"   "16:670964"   "16:783865"
  [6] "16:1129384"  "16:1536499"  "16:1634385"  "16:1890373"  "16:1918125"
 [11] "16:2018650"  "16:2026986"  "16:3640784"  "16:3656625"  "16:4487486"
 [16] "16:7772926"  "16:9382304"  "16:11541896" "16:11876203" "16:15131974"
 [21] "16:15962510" "16:16276292" "16:16278869" "16:18961444" "16:20435314"
 [26] "16:20446192" "16:20791506" "16:20966362" "16:20975505" "16:20976360"
 [31] "16:22038566" "16:22059169" "16:23227396" "16:23811872" "16:24801468"
 [36] "16:24802325" "16:28865042" "16:28867804" "16:29830829" "16:29835588"
 [41] "16:29857316" "16:29860685" "16:30010475" "16:30094180" "16:34864943"
 [46] "16:48177220" "16:49671218" "16:52738113" "16:56861444" "16:56930251"
 [51] "16:57101373" "16:57474687" "16:57927121" "16:57927336" "16:57928156"
 [56] "16:57953715" "16:57954986" "16:57955613" "16:57955692" "16:57956231"
 [61] "16:58061026" "16:58064274" "16:58072153" "16:58093750" "16:58095031"
 [66] "16:58096591" "16:58100228" "16:58101230" "16:58314598" "16:58629135"
 [71] "16:62369902" "16:63283562" "16:66586189" "16:66589440" "16:66593067"
 [76] "16:66604604" "16:67114330" "16:67304915" "16:67320223" "16:67409180"
 [81] "16:68261092" "16:70978985" "16:74750351" "16:75261647" "16:75262639"
 [86] "16:75266918" "16:75283093" "16:75288237" "16:75288674" "16:75301196"
 [91] "16:75301232" "16:75388731" "16:75434558" "16:75438777" "16:75493481"
 [96] "16:75534058" "16:75536272" "16:77930075" "16:81128230" "16:81161569"
[101] "16:81187709" "16:83246883" "16:88116843" "16:88973220"
> Data1[1:10,1:10]
   X1.10759741_T X2.63375893_G X1.82944926_C X2.142291951_C X1.89628674_G
1              1             0             0              1             0
2              0             0             1              0             0
3              0             0             0              2             0
4              0             0             0              0             0
5              0             0             2              0             1
6              0             0             1              2             0
7              0             0             0              1             0
8              0             0             1              1             0
9              0             0             0              0             0
10             0             0             0              1             0
   X1.238625898_T X1.30113436_G X2.70564002_A X1.171543747_C X1.192618028_C
1               0             1             1              0              0
2               1             1             0              0              1
3               1             0             2              0              0
4               1             0             0              0              0
5               0             1             1              0              0
6               0             0             1              0              0
7               2             0             0              0              0
8               1             2             0              0              0
9               1             0             1              0              0
10              0             0             1              0              0
> Data2[1:10,1:10]
      X1.10759741_T X2.63375893_G X1.82944926_C X2.142291951_C X1.89628674_G
 [1,]             0             0             0              1             0
 [2,]             0             0             1              0             0
 [3,]             0             0             0              0             1
 [4,]             0             0             0              1             0
 [5,]             1             0             2              1             0
 [6,]             0             0             0              1             0
 [7,]             0             0             0              0             0
 [8,]             0             0             0              0             0
 [9,]             0             0             0              1             0
[10,]             0             0             0              1             0
      X1.238625898_T X1.30113436_G X2.70564002_A X1.171543747_C X1.192618028_C
 [1,]              0             1             1              0              0
 [2,]              1             0             1              0              0
 [3,]              1             0             2              0              0
 [4,]              2             0             0              0              0
 [5,]              1             0             0              0              0
 [6,]              0             1             0              0              0
 [7,]              1             2             0              0              0
 [8,]              0             0             0              0              0
 [9,]              2             1             1              0              0
[10,]              2             1             1              0              0



~~~

