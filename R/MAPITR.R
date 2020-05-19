#' @title MArginal ePIstasis Test for Regions
#' (\code{MAPITR})
#'
#' @description Run \code{MAPITR} on a set of phenotypes that each have
#' univariate GWAS statistics on the same set of SNPs
#' 
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes. No default value.
#' 
#' @param ... Additional optional arguments.
#' 
#' @return A list containing model, SNP, and posterior information for
#' both the previously significant univariate SNPs (\code{PreviousSNPs})
#' and the newly significant multivariate SNPs (\code{NewSNPs}). For a 
#' full breakdown of the \code{bmass} output list structure, please see 
#' the associated vignettes.
#'
#' @examples
#' Phenotypes <- c("bmass_SimulatedData1", "bmass_SimulatedData2")
#' bmassOutput <- bmass(Phenotypes, bmass_SimulatedSigSNPs)
#' summary(bmassOutput)
#' bmassOutput$NewSNPs$SNPs
#'
#' @section Other Examples:
#' \code{bmass(c("HDL","LDL","TG","TC"), GWASsnps, NminThreshold = 50000)}
#' \code{bmass(c("HDL","LDL","TG","TC"), GWASsnps, GWASThreshValue = 1e-8,
#'   NminThreshold = 50000, PrintProgress = TRUE)} 
#' \code{bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, GWASThreshFlag = FALSE,
#'   SNPMarginalUnivariateThreshold = 1e-4,
#'   SNPMarginalMultivariateThreshold = 1e-4,
#'   PrintMergedData = TRUE)} 
#' \code{bmassOutput <- bmass(c("HDL","LDL","TG","TC"),
#'   GWASsnps, NminThreshold = 50000)} 
#'
#' @export
#' 
MAPITR <- function (Phenotypes, Genotypes, Pathways, Covariates = NULL, CenterStandardize = TRUE, RegressPhenotypes = TRUE)
	return(MAPITRmain(Phenotypes, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, ...))
}

#DataSources, GWASsnps = NULL, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, GWASThreshFlag = TRUE, GWASThreshValue = 5e-8, NminThreshold = 0, PrintMergedData = FALSE, PrintProgress = FALSE, ...) {
#DataSources, GWASsnps, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, GWASThreshFlag, GWASThreshValue, NminThreshold, PrintMergedData, PrintProgress, MergedDataSources = NULL, ZScoresCorMatrix = NULL, ExpectedColumnNames = c("Chr", "BP", "Marker", "MAF", "A1", "Direction", "pValue", "N"), GWASsnps_AnnotateWindow = 5e5, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), ProvidedPriors = NULL, UseFlatPriors = NULL, PruneMarginalSNPs = TRUE, PruneMarginalSNPs_bpWindow = 5e5, bmassSeedValue = 1) {

MAPITRmain <- function (PhenotypesVector, Genotypes, Pathways, GRM_Grand = Null, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) {

        MAPITRprocessing <- list()
	MAPITRoutput <- list()
	MAPITRoutput$LogFile <- c();
	MAPITRoutput$pValue <- NULL
	MAPITRoutput$PVE <- NULL
	cores = detectCores()

        MAPITRoutput$LogFile <- rbind(MAPITRoutput$LogFile, paste(format(Sys.time()), " -- beginning MAPITR.", sep=""))
        if (PrintProgress == TRUE) {
                write(paste(format(Sys.time()), " -- beginning MAPITR.", sep=""), stderr())
        }

	#DataChecks
	MAPITRoutput$LogFile <- DataChecks(PhenotypesVector, Genotypes, Pathways, Covariates, MAPITRoutput$LogFile)
#	MAPITRprocessing$log <- DataChecks(PhenotypesVector, Genotypes, Pathways, Covariates

#	                bmassOutput$LogFile <- CheckIndividualDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
#	MAPITRoutput[c("MergedDataSources", "LogFile")] <- MergeDataSources(DataSources, MAPITRoutput$LogFile)[c("MergedDataSources", "LogFile")]


	MAPITRoutput.temp1 <- PreprocessData(PhenotypesVector, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, MAPITRoutput$LogFile)
	PhenotypesMatrix <- MAPITRoutput.temp1$PhenotypesMatrix
	Genotypes <- MAPITRoutput.temp1$Genotypes
	Pathways.Full <- MAPITRoutput.temp1$Pathways.Full
	MAPITRoutput$LogFile <- MAPITRoutput.temp1$LogFile
	rm(MAPITRoutput.temp1)	

	if ( ) {
		MAPITRoutput.temp2 <- RunMAPITR.NothingProvided(PhenotypesVector, Genotypes, Pathways.Full, cores, MAPITRoutput$LogFile) 
	} elsif ( ) { 

	} else if
		

	return(MAPITRoutput) 

}


