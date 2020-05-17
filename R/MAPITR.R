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
MAPITROld <- function (DataSources, GWASsnps = NULL, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, GWASThreshFlag = TRUE, GWASThreshValue = 5e-8, NminThreshold = 0, PrintMergedData = FALSE, PrintProgress = FALSE, ...) {
	return(bmassMain(DataSources, GWASsnps, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, GWASThreshFlag, GWASThreshValue, NminThreshold, PrintMergedData, PrintProgress, ...))
}

MAPITR <- function (DataSources, GWASsnps, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, GWASThreshFlag, GWASThreshValue, NminThreshold, PrintMergedData, PrintProgress, MergedDataSources = NULL, ZScoresCorMatrix = NULL, ExpectedColumnNames = c("Chr", "BP", "Marker", "MAF", "A1", "Direction", "pValue", "N"), GWASsnps_AnnotateWindow = 5e5, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), ProvidedPriors = NULL, UseFlatPriors = NULL, PruneMarginalSNPs = TRUE, PruneMarginalSNPs_bpWindow = 5e5, bmassSeedValue = 1) {

        bmassOutput <- list()
	bmassOutput$MergedDataSources <- NULL
        if (!is.null(MergedDataSources)) {
		bmassOutput$MergedDataSources <- MergedDataSources
	}
	bmassOutput$ZScoresCorMatrix <- NULL
        if (!is.null(ZScoresCorMatrix)) {
		bmassOutput$ZScoresCorMatrix <- ZScoresCorMatrix
	}
	bmassOutput$MarginalSNPs <- list()
	bmassOutput$Models <- NULL
	bmassOutput$ModelPriors <- NULL
        bmassOutput$PreviousSNPs <- list()
        bmassOutput$NewSNPs <- list()
        bmassOutput$GWASlogBFMinThreshold <- NULL
	bmassOutput$LogFile <- c()
	bmassOutput$MarginalSNPs <- list()
	bmassOutput$Models <- NULL
	bmassOutput$ModelPriors <- NULL
        bmassOutput$PreviousSNPs <- list()
        bmassOutput$NewSNPs <- list()
        bmassOutput$GWASlogBFMinThreshold <- NULL
	bmassOutput$LogFile <- c()

        bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- beginning bmass.", sep=""))
	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- beginning bmass.", sep=""), stderr())
	}

        #Loading and checking data
        #~~~~~~

	if (!is.null(MergedDataSources)) {
		bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- MergedDataSources was provided, skipping merging data step.", sep=""))
		if (PrintProgress == TRUE) {
        		write(paste(format(Sys.time()), " -- MergedDataSources was provided, skipping merging data step.", sep=""), stderr())
		}
	} else {
		if (PrintProgress == TRUE) {
        		write(paste(format(Sys.time()), " -- Checking individual datasource files and merging datasets.", sep=""), stderr())
		}
		bmassOutput$LogFile <- CheckIndividualDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
		
		bmassOutput[c("MergedDataSources", "LogFile")] <- MergeDataSources(DataSources, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]
	}

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Annotating merged datasources with GWAS SNPs (if provided) then processing finalized merged dataset.", sep=""), stderr())
	}
	bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]

	bmassOutput[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")] <- ProcessMergedAndAnnotatedDataSources(DataSources, bmassOutput$MergedDataSources, bmassOutput$ZScoresCorMatrix, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, bmassOutput$LogFile)[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")]

	if (PrintMergedData == FALSE) {
		bmassOutput$MergedDataSources <- NULL
	}

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Calculating logBFs and posteriors (if applicable).", sep=""), stderr())
	}
	
	bmassOutput[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] <- GetLogBFsFromData(DataSources, bmassOutput$MarginalSNPs, bmassOutput$ZScoresCorMatrix, SigmaAlphas, bmassOutput$LogFile)[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] 

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")] <- DetermineAndApplyPriors(DataSources, bmassOutput$MarginalSNPs, GWASsnps, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")]
	
	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Getting final list of MarginalSNPs, PreviousSNPs, and NewSNPs (where applicable).", sep=""), stderr())
	}

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")] <- FinalizeAndFormatResults(DataSources, bmassOutput$MarginalSNPs, bmassOutput$PreviousSNPs, GWASsnps, bmassOutput$GWASlogBFMinThreshold, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, NminThreshold, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")]

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Finishing bmass, exiting.", sep=""), stderr())
	}

	return(bmassOutput)
}
