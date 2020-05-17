CheckNumericClass <- function (InputData) {
	returnValue <- FALSE
	
	if (is.numeric(InputData)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

CheckForNAs <- function (InputData) {
	returnValue <- FALSE
	
	if (is.na(InputData)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

DataChecks <- function (PhenotypesVector, Genotypes, Pathways, Covariates, LogFile) {

	#Checking for NAs in phenotype vector
	#unit test
	if (is.na(PhenotypesVector)) {
		stop(Sys.time(), " -- there are NAs in the phenotype vector. Please remove and align the remaining files (eg genotype matrix).");
	}
	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- PhenotypesVector passed NA check.", sep=""))

	#Checking for missingness in genotype matrix
	#unit test
	if (apply(Genotypes, 2, is.na)) {
		stop(Sys.time(), " -- there are NAs in the genotype matrix. There must be zero missingness in the genotype matrix. Please correct and rerun.");
	}
	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Genotypes passed NA check.", sep=""))

#	#Checking for NAs in pathways file
#	if (apply(Pathways, 1, 

#	#If using covariates, checking for NAs in covariate matrix

	return(LogFile)

}

PreprocessData <- function (PhenotypesVector, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, LogFile) {

	PreprocessData.Output <- list();
	PhenotypesMatrix <- c();

	#unit test this
	if (CenterStandardize == TRUE) {
		Genotypes.Mean <- apply(Genotypes, 2, mean); 
		Genotypes.SD <- apply(Genotypes, 2, sd); 
		Genotypes <- t((t(Genotypes)-Genotypes.Mean)/Genotypes.SD);
	}
	Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \i

	#unit test??
	if (RegressPhenotypes == TRUE) { 
		for (i in 1:nrow(Pathways)) { 
			Genotypes.Pathway <- Genotypes[,Pathway];
			PhenotypesMatrix <- cbind(PhenotypesMatrix, residuals(lm(PhenotypesVector ~ Genotypes.Pathway - 1)))	
		}
	} else {
		for (i in 1:nrow(Pathways)) { 
			PhenotypesMatrix <- cbind(PhenotypesMatrix, PhenotypesVector);
		}
	}

	return(list(PhenotypesMatrix=PhenotypesMatrix, Genotypes=Genotypes, LogFile=LogFile))
	return(PreprocessData.Output);
	return(list(MergedDataSources=MergedDataSources, MarginalSNPs=MarginalSNPs, ZScoresCorMatrix=ZScoresCorMatrix, LogFile=LogFile))

}
