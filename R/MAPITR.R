#' @title MArginal ePIstasis Test for Regions
#' (\code{MAPITR})
#'
#' @description Run \code{MAPITR} for a group of pathways on a single
#' phenotype and a set of genome-wide SNPs
#'
#' @param Phenotype
#'
#' @param Genotypes
#'
#' @param Pathways
#'
#' @param Covariates
#'
#' @param CenterStandardize
#'
#' @param RegressPhenotypes
#'
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes. No default value.
#' 
#' @param ... Additional optional arguments.
#'
#' @return A matrix containing in the first column the list of pathways
#' that were analyzed, in the second column the associated 
#' \code{MAPITR} p-values for each pathway, and in the third column the
#' associated \code{MAPITR} PVEs for each pathway.
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
MAPITR <- function (Phenotype, Genotypes, Pathways, Covariates = NULL, CenterStandardize = TRUE, RegressPhenotypes = TRUE)
	return(MAPITRmain(Phenotype, Genotypes, Pathways, Covariates, CenterStandardize, RegressPhenotypes, ...))
}

MAPITRmain <- function (Phenotype, Genotypes, Pathways, GRM_Grand = NULL, GRM_Pathway = NULL, Covariates, CenterStandardize = TRUE, RegressPhenotypes = TRUE, PrintProgress = FALSE) {

        MAPITRprocessing <- list()
	MAPITRoutput <- list()
	MAPITRoutput$LogFile <- c()
	MAPITRoutput$pValues <- NULL
	MAPITRoutput$PVE <- NULL
	cores = detectCores()

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

	#Running appropriate version of MAPITR
	if (is.null(Covariates)) {
		MAPITRoutput.temp2 <- RunMAPITR.Base(PhenotypeMatrix, Genotypes, Pathways.Full, cores, MAPITRoutput$LogFile) 
		MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp2$Est, MAPITRoutput.temp2$Eigenvalues)
		rm(MAPITRoutput.temp2)
	} else if (!is.null(Covariates)) { 
		MAPITRoutput.temp3 <- RunMAPITR.wCovs(PhenotypeMatrix, Genotypes, Covariates, Pathways.Full, cores, MAPITRoutput$LogFile) 
		MAPITRoutput$pValues <- GetMAPITRpValues(MAPITRoutput.temp3$Est, MAPITRoutput.temp3$Eigenvalues)
		rm(MAPITRoutput.temp3)
	} else {
		stop(Sys.time(), " -- 'Covariates' is neither null or not null. This should not happen. Contact current maintainer of code.")
	}
		
	#Postprocessing of results (if needed)
	Pathway.Names <- Pathways[,1] 
	MAPITRoutput.Final <- cbind(Pathway.Names, MAPITRoutput$pValues, MAPITRoutput$PVE) 

	return(MAPITRoutput.Final) 

}
