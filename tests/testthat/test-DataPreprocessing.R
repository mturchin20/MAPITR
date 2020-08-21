context("Tests for DataPreprocessing.R")

assign("Pathways_Example", c(1,2,3,4), envir = .GlobalEnv)
assign("CenterStandardize_Orig", matrix(c(1,2,0,1,0,0,2,1,1,1,2,1), ncol=3, byrow=FALSE), envir = .GlobalEnv)
assign("CenterStandardize_Mean", c(1.00,0.75,1.25), envir = .GlobalEnv)
assign("CenterStandardize_SD", c(0.8164966,0.9574271,0.5000000), envir = .GlobalEnv)
assign("CenterStandardize_Final", matrix(c(0.0000000,1.2247449,-1.2247449,0.0000000,-0.7833495,-0.7833495,1.3055824,0.2611165,-0.5000000,-0.5000000,1.5000000,-0.5000000), ncol=3, byrow=FALSE), envir = .GlobalEnv)
assign("Residuals_Example",
envir = .GlobalEnv)



Pathways.Full <- lapply(strsplit(as.character(Pathways[,2]), ","), as.numeric)
Genotypes.Mean <- apply(Genotypes, 2, mean);
                Genotypes.SD <- apply(Genotypes, 2, sd);
                Genotypes <- t((t(Genotypes)-Genotypes.Mean)/Genotypes.SD);
 Genotypes.Pathway <- Genotypes[,as.numeric(unlist(Pathways.Full[i]))];
                        PhenotypeMatrix <- cbind(PhenotypeMatrix, residuals(lm(as.matrix(PhenotypesVector) ~ as.matrix(Genotypes.Pathway) - 1)))

test_that("Checking behavior of pathway format splitting and setup", {
	 expect_equal(lapply(strsplit(as.character("1,2,3,4"), ","), as.numeric)[[1]], Pathways_Example)
})

test_that("Checking behavior of centering and standardizing genotypes procedure", {
	expect_equal(apply(CenterStandardize_Orig, 2, mean), CenterStandardize_Mean)
	expect_equal(apply(CenterStandardize_Orig, 2, sd), CenterStandardize_SD, tolerance=1e-6)
	expect_equal(t((t(CenterStandardize_Orig) - CenterStandardize_Mean)/CenterStandardize_SD), CenterStandardize_Final, tolerance=1e-6)
}

test_that("Checking behavior of linear regression residuals extraction", {

}


test_that("CollapseSigmaAlphasTogether sums multiple entries of the same 'model' over different sigma_alpha hyperparameter values", {
	expect_equal(CollapseSigmaAlphasTogether(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), 3), CollapseTest1)
})

test_that("CheckForAndReplaceOnes checks for 1 entries and replaces them with 0s", {
	expect_equal(c(apply(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), c(1,2), CheckForAndReplaceOnes)), CheckOnesTest1)
	expect_equal(c(apply(10^matrix(0, ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)), CheckOnesTest2) 
})

test_that("CheckForAndReplaceZeroes checks for 0 entries and replaces them with 1s", {
	expect_equal(c(apply(apply(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), c(1,2), CheckForAndReplaceOnes), c(1,2), CheckForAndReplaceZeroes)), CheckZeroesTest1)
	expect_equal(c(apply(10^matrix(c(0,1), ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)), CheckZeroesTest2) 
})

test_that("GetSumAcrossSigmaAlphas_withPriors sums multiple entries (in log10) of the same 'model' over different sigma_alpha hyperparameter values while multiplying each entry by that model's prior", {
	expect_equal(c(GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(MarginalSNPs_logBFs_Stacked)), nrow=length(ModelPriors), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), 3, 3)), CollapseTest_wPriors, tolerance=1e-6)
	expect_equal(c(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)[seq.int(1, by=2, length.out=3),]), CollapseTest_wPriors_Basic1) 
})

rm(CollapseTest1, CheckOnesTest1, CheckOnesTest2, CheckZeroesTest1, CheckZeroesTest2, ModelPriors, MarginalSNPs_logBFs_Stacked, CollapseTest_wPriors, CollapseTest_wPriors_Basic1, envir = .GlobalEnv)



rm(
, envir = .GlobalEnv)
