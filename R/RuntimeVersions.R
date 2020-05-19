#' @useDynLib MAPITR
#' @import doParallel
#' @import Rcpp
#' @import RcppArmadillo
#' @import RcppParallel
RunMAPITR.Base <- function (PhenotypesMatrix, Genotypes, Pathways.Full, cores, LogFile) {

	RunMAPITR.Base.Output <- list()

	#MAPITR expects a n x r phenotype matrix, a p x n genotype matrix, and a list of SNP indices for each pathway
	RunMAPITR.Base.Output.temp1 <- MAPITR(as.matrix(PhenotypesMatrix),t(as.matrix(Genotypes)),Pathways.Full,cores=cores)

	return(list(Est=RunMAPITR.Base.Output.temp1$Est, Eigenvalues=RunMAPITR.Base.Output.temp1$Eigenvalues, PVE=RunMAPITR.Base.Output.temp1$PVE, LogFile=LogFile))

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

#' @importFrom stats pnorm runif
DetermineAndApplyPriors <- function(DataSources, MarginalSNPs, GWASsnps, SigmaAlphas, Models, ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, LogFile) {

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, ModelPriors=ModelPriors_Used, GWASlogBFMinThreshold=PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, LogFile=LogFile))

}
