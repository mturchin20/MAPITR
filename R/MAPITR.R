#' @title MArginal ePIstasis Test for Regions
#' (\code{MAPITR})
#'
#' @description Run \code{MAPITR} for a group of pathways on a single
#' phenotype and a set of genome-wide SNPs
#'
#' @param Phenotype A vector containing phenotypic values for all
#' individuals being analyzed. No default value.
#'
#' @param Genotypes A n x p matrix containing the genotypes (0/1/2) for
#' all p SNPs across all n individuals. No default value.
#'
#' @param Pathways A r x 2 matrix containing the pathway names and then
#' a comma-separated list of each \code{Genotypes} column index
#' representing each SNP in the associated pathway. Note, this second 
#' column of comma-separated indices are the numeric positions for each 
#' SNP in \code{Genotypes} and not the SNP IDs or column names. No
#' default value.
#'
#' @param Covariates A n x q matrix containing any q additional covariates
#' that should be included in the M-projection matrix of the model. 
#' See Turchin et al. 2019 for details. Note that these are covariates
#' which are applied to both sides of the model, ie the phenotype as well
#' as the genotypes. A y-intercept term is automatically included and does
#' not need to be part of this n x q matrix. No default value.
#'
#' @param CenterStandardize A logical \code{TRUE}/\code{FALSE} flag that
#' indicates whether the genotype matrix \code{Genotypes} should be
#' centered and standardized before analysis. This is a recommended step.
#' Indicate \code{FALSE} if this is a preprocessing step that has
#' already been done prior to running \code{MAPITR}. The default value 
#' is TRUE.  
#'
#' @param ... Additional optional arguments.
#'
#' @return A matrix containing in the first column the list of pathways
#' that were analyzed, in the second column the associated 
#' \code{MAPITR} p-values for each pathway, and in the third column the
#' associated \code{MAPITR} PVEs for each pathway.
#' 
#' @examples
#' Phenotype <- c("MAPITR_SimulatedPhenotype")
#' Genotypes <- c("MAPITR_SimulatedGenotypes")
#' Pathways <- c("MAPITR_SimulatedPathways")
#' MAPITROutput <- bmass(MAPITR_SimulatedPhenotype, MAPITR_SimulatedGenotypes, MAPITR_SimulatedPathways)
#' summary(MAPITROutput)
#' MAPITR[MAPITR$pValues < 1e-4,1]
#'
#' @export
#' 
MAPITR <- function (Phenotype, Genotypes, Pathways, Covariates = NULL, CenterStandardize = TRUE) {
	return(MAPITRmain(Phenotype, Genotypes, Pathways, Covariates, CenterStandardize, ...))
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
