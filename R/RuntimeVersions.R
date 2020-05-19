#' @useDynLib MAPITR
#' @import doParallel
#' @import Rcpp
#' @import RcppArmadillo
#' @import RcppParallel
RunMAPITR.NothingProvided <- function (PhenotypesMatrix, Genotypes, Pathways.Full, cores, LogFile) {

	RunMAPITR.NothingProvided.Output <- list()

	#MAPITR expects a n x r phenotype matrix, a p x n genotype matrix, and a list of SNP indices for each pathway
	RunMAPITR.NothingProvided.Output.temp1 <- MAPITR(as.matrix(PhenotypesMatrix),t(as.matrix(Genotypes)),Pathways.Full,cores=cores)

	return(list(Est=RunMAPITR.NothingProvided.Output.temp1$Est, Eigenvalues=RunMAPITR.NothingProvided.Output.temp1$Eigenvalues, PVE=RunMAPITR.NothingProvided.Output.temp1$PVE, LogFile=LogFile))

}


RunMAPITR.wCovs <- function (Phenotypes, Genotypes, Pathway, Covariates, CenterStandardize) {

	RunMAPITR.wCovs.Output <- list()
	
	#MAPITR.wCovs expects a n x r phenotype matrix, a p x n genotype matrix, a z x n covariate matrix, and a list of SNP indices for each pathway
	RunMAPITR.wCovs.Output.temp1 <- MAPITR.wCovs(as.matrix(PhenotypesMatrix),t(as.matrix(Genotypes.Pathway)),Pathways.Full,t(as.matrix(Covariates)),cores=cores)

	return(RunMAPITR.wCovs.Output)

}

#RunMAPITR.GRMsProvided <- function ( ) {
#
#	RunMAPITR.GRMsProvided.Output <- list()
#
#	return(RunMAPITR.GRMsProvided.Output)
#
#}

#' @importFrom CompQuadForm davies
GetMAPITRpValues <- function (Est, Eigenvalues, acc=1e-8) {

	pValues <- c()
	for (i in 1:length(Est)) { 
		Lambda <- sort(Eigenvalues[,i], decreasing=TRUE)
		Davies.Output <- davies(Est[i], lambda=Lambda, acc=acc)
		pValues <- c(pValues, 2*min(1-Davies.Output$Qq, Davies.Output$Qq))
	}

	return(pValues)

}

#Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
#Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
#pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \




# This function expects columns from nSigmaAlphas stacked Model x SNP
# matrices of posterior probabilities, such that a single column (ie
# single SNP) is converted to a Model x nSigmaAlphas matrix and summed
# across rows for a single posterior probability per Model (eg the
# "marginal" of across all sigma_alphas)
CollapseSigmaAlphasTogether <- function (inputValues1, nSigmaAlphas) {
        CollapsedInputs <- apply(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE), 1, sum)
        return(CollapsedInputs)
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

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, ModelPriors=ModelPriors_Used, GWASlogBFMinThreshold=PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, LogFile=LogFile))

}
