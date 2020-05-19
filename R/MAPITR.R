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

MAPITRmain <- function (PhenotypesVector, Genotypes, Pathways, GRM_Grand = Null, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) {

        MAPITRprocessing <- list()
	MAPITRoutput <- list()
	MAPITRoutput$LogFile <- c();
	MAPITRoutput$pValues <- NULL
	MAPITRoutput$PVE <- NULL
	cores = detectCores()

        MAPITRoutput$LogFile <- rbind(MAPITRoutput$LogFile, paste(format(Sys.time()), " -- beginning MAPITR.", sep=""))
        if (PrintProgress == TRUE) {
                write(paste(format(Sys.time()), " -- beginning MAPITR.", sep=""), stderr())
        }

	#Data Checks
	MAPITRoutput$LogFile <- DataChecks(PhenotypesVector, Genotypes, Pathways, Covariates, MAPITRoutput$LogFile)

	#Preprocessing Data
	MAPITRoutput.temp1 <- PreprocessData(PhenotypesVector, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, MAPITRoutput$LogFile)
	PhenotypesMatrix <- MAPITRoutput.temp1$PhenotypesMatrix
	Genotypes <- MAPITRoutput.temp1$Genotypes
	Pathways.Full <- MAPITRoutput.temp1$Pathways.Full
	MAPITRoutput$LogFile <- MAPITRoutput.temp1$LogFile
	rm(MAPITRoutput.temp1)	

	#Running appropriate version of MAPITR
	
	if (is.null(Covariates)) {
		MAPITRoutput.temp2 <- RunMAPITR.Base(PhenotypesVector, Genotypes, Pathways.Full, cores, MAPITRoutput$LogFile) 
		MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp2$Est, MAPITRoutput.temp2$Eigenvalues)
	} else if (!is.null(Covariates)) { 
		MAPITRoutput.temp3 <- RunMAPITR.wCovs(PhenotypesVector, Genotypes, Covariates, Pathways.Full, cores, MAPITRoutput$LogFile) 
		MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp3$Est, MAPITRoutput.temp3$Eigenvalues)
	} else {
		stop(Sys.time(), " -- there are NAs in the genotype matrix. There must be zero missingness in the genotype matrix. Please correct and rerun.");	
	}
		
	#Postprocessing of results (if needed)
	Pathway.Names <- Pathways[,1]; 
	MAPITRoutput.Final <- cbind(Pathway.Names, MAPITRoutput$pValues, MAPITRoutput$PVE); 

	return(MAPITRoutput.Final) 

}
