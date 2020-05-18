#' @useDynLib MAPITR
#' @import doParallel
#' @import Rcpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @import CompQuadForm
RunMAPITR.NothingProvided <- function (PhenotypesMatrix, Genotypes, Pathways.Full, LogFile) {

	RunMAPITR.NothingProvided.Output <- list()

	#MAPITR expects a n x r phenotype matrix, a p x n genotype matrix, and a list of SNP indices for each pathway
	RunMAPITR.NothingProvided.Output.temp1 <- MAPITR(PhenotypesMatrix,t(Genotypes.Pathway),Pathways.Full,cores=cores)
#	List MAPITR(mat Y,mat X,mat regions,int cores = 1){
#	InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,t(as.matrix(Z)),Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores);

	return(list(Est=RunMAPITR.NothingProvided.Output.temp1$Est, Eigenvalues=RunMAPITR.NothingProvided.Output.temp1$Eigenvalues, PVE=RunMAPITR.NothingProvided.Output.temp1$PVE, LogFile=LogFile))
#	return(list(PhenotypesMatrix=PhenotypesMatrix, Genotypes=Genotypes, LogFile=LogFile))

}

sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp"); 

 Y.Check.Pheno.noNAs <- Y.Check.Pheno[neg.is.na(Y.Check.Pheno)];
 
Pathways.Regions <- list(); 
cores = detectCores(); 
InterPath.output <- list(); 
InterPath.output$Est <- c();
 InterPath.output$Eigenvalues <- c(); 
InterPath.output$PVE <- c(); \

Pathways.Regions[[1]] 
                                       
K <- 1/ncol(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \

InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,t(as.matrix(Z)),Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); 

InterPath.output$Est <- c(InterPath.output$Est, InterPath.output.temp$Est); InterPath.output$Eigenvalues <- cbind(InterPath.output$Eigenvalues, InterPath.output.temp$Eigenvalues); InterPath.output$PVE <- c(InterPath.output$PVE, InterPath.output.temp$PVE); 
                                
Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \


	return(RunMAPITR.NothingProvided.Output)

}

RunMAPITR.wCovs <- function (Phenotypes, Genotypes, Pathway, Covariates, CenterStandardize) {

	RunMAPITR.wCovs.Output <- list()
	
	RunMAPITR.NothingProvided.Output.temp1 <- MAPITR(PhenotypesMatrix,t(Genotypes.Pathway),PhenotypesMatrix,as.matrix(GRM_Grand),as.matrix(GRM_Patht(as.matrix(Z)),cores=cores)

	return(RunMAPITR.wCovs.Output)

}

RunMAPITR.GRMsProvided <- function ( ) {

	RunMAPITR.GRMsProvided.Output <- list()

	return(RunMAPITR.GRMsProvided.Output)

}




# This function expects columns from nSigmaAlphas stacked Model x SNP
# matrices of posterior probabilities, such that a single column (ie
# single SNP) is converted to a Model x nSigmaAlphas matrix and summed
# across rows for a single posterior probability per Model (eg the
# "marginal" of across all sigma_alphas)
CollapseSigmaAlphasTogether <- function (inputValues1, nSigmaAlphas) {
        CollapsedInputs <- apply(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE), 1, sum)
        return(CollapsedInputs)
}

CheckForAndReplaceZeroes <- function(x) {
        returnValue1 <- x
        if (x == 0) {
                returnValue1 <- 1
        }
        return(returnValue1)
}

#This function is expecting nSigmaAlphas Model x SNP matrices of
#logBFs stacked ontop of one another. ModelPriors_Matrix is a matrix
#containing the vector of model priors (ModelPriors = one set of model
#priors * nSigmaAlphas) replicated as a column for each snp.
GetSumAcrossSigmaAlphas_withPriors <- function(logBFs1, ModelPriors_Matrix, nGammas, nSigmaAlphas) {
	WeightedSumAcrossAlphaSigmas <- matrix(0, ncol=ncol(logBFs1), nrow=nGammas)
        for (i in 1:nGammas) {
                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
                max <- apply(logBFs1[SigmaAlpha_Coordinates,], 2, max)
                logBFs1[SigmaAlpha_Coordinates,] <- logBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(logBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(logBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
		WeightedSumAcrossAlphaSigmas[i,] <- log10(sapply(apply(ModelPriors_Matrix[SigmaAlpha_Coordinates,] * apply(10^logBFs1[SigmaAlpha_Coordinates,], c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + max
        }
        return(WeightedSumAcrossAlphaSigmas)
}

#' @importFrom stats pnorm runif
DetermineAndApplyPriors <- function(DataSources, MarginalSNPs, GWASsnps, SigmaAlphas, Models, ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, LogFile) {

	MarginalSNPs_logBFs_Stacked <- MarginalSNPs$logBFs	
	MarginalSNPs_logBFs_Stacked_AvgwPrior <- NULL
        ModelPriors_Used <- ModelPriors
		
	PreviousSNPs <- list()
        PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- NULL
	ZScoreHitFlag1 <- c()
	if (GWASThreshFlag) {
		ZScoreHitFlag1 <- rep(0, nrow(MarginalSNPs$SNPs))
		ZScoreHitFlag1[2*pnorm(apply(abs(MarginalSNPs$SNPs[,grep("ZScore", colnames(MarginalSNPs$SNPs))]),1,max),0,1,lower.tail=FALSE) < GWASThreshValue] <- 1
	}

	if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, ProvidedPriors)
                ModelPriors_Used <- ProvidedPriors
		PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked
        } else if (is.null(GWASsnps) && UseFlatPriors == TRUE) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))

		Prior_FlatUnif <- normalize(rep(c(0,ModelPriors[-1]),length(SigmaAlphas)))

                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_FlatUnif)
                ModelPriors_Used <- Prior_FlatUnif

        } else if (is.null(GWASsnps) && UseFlatPriors == FALSE) {
		stop("Error 1a (GetResultsFromMarginalSNPsAndFormat.R) -- code does not currently accept 'UseFlatPrior == FALSE' while not providing a list of GWAS snps (ie is.null(GWASsnps) == TRUE); please rerun with different argument options"); 
	} else {
	
		if (!is.null(bmassSeedValue)) {
			set.seed(bmassSeedValue)
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- setting seed with the following value: ", bmassSeedValue, ".", sep=""))
		}
		else {
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- no seed value via bmassSeedValue provided.", sep=""))
		}
		
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
               
		PreviousSNPs_logBFs_Stacked <- c() #Will be matrix of nSigmaAlphas x nSNPs
		if (GWASThreshFlag) {
			PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==1]) 
		} 
		else if (!GWASThreshFlag) {
			PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1]) 
		} 
		else {
			stop("Error 2a (GetResultsFromMarginalSNPsAndFormat.R) -- Unexpected conditional outcome, 'if (GWASThreshFlag)' nor 'if (!GWASThreshFlag)' are true");
		}

                Prior_PreviousSNPsEB <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors, 100) #Vector with nModels*nSigmaAlphas entries
                Prior_PreviousSNPsEB_check2 <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors*runif(length(ModelPriors)), 100)

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_PreviousSNPsEB)
                ModelPriors_Used <- Prior_PreviousSNPsEB
                
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked
        }

        if (is.null(MarginalSNPs_logBFs_Stacked_AvgwPrior)) {
                stop(Sys.time(), " -- No average log BFs were returned from method. Check if all input variables are as the method expects.")
        }
        
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked
        MarginalSNPs$SNPs$logBFWeightedAvg <- MarginalSNPs_logBFs_Stacked_AvgwPrior
		
	if (!is.null(GWASsnps)) {
		if (GWASThreshFlag) {
			PreviousSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==1,]
			PreviousSNPs$DontPassSNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==0,]
		} else if (!GWASThreshFlag) {
			PreviousSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1,]
		} else {
			stop("Error 2b (GetResultsFromMarginalSNPsAndFormat.R) -- Unexpected conditional outcome, 'if (GWASThreshFlag)' nor 'if (!GWASThreshFlag)' are true");
		}	
	
		if (dim(PreviousSNPs$SNPs)[1] > 0) {
			PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- min(PreviousSNPs$SNPs$logBFWeightedAvg)
		}
	}

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, ModelPriors=ModelPriors_Used, GWASlogBFMinThreshold=PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, LogFile=LogFile))

}
