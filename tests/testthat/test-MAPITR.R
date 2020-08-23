context("Tests for MAPITR.R") 

data(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways)

assign("MAPITR_TestData_Genotypes", MAPITR_TestData_Genotypes, envir = .GlobalEnv)
assign("MAPITR_TestData_Phenotype", MAPITR_TestData_Phenotype, envir = .GlobalEnv)
assign("MAPITR_TestData_Pathways", MAPITR_TestData_Pathways, envir = .GlobalEnv)
assign("MAPITR_Output", MAPITR(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways), envir = .GlobalEnv)
assign("MAPITR_Output_pValue1", c(0.019334), envir = .GlobalEnv)
assign("MAPITR_Output_Est1", c(1.07617), envir = .GlobalEnv)
assign("MAPITR_Output_PVE1", c(1.2418), envir = .GlobalEnv)
assign("MAPITR_Output_Eigen1", c(-0.01903), envir = .GlobalEnv)
assign("MAPITR_Output_Eigen2", c(-0.01732), envir = .GlobalEnv)

test_that("MAPITR runs the main MAPITR function", {
	expect_equal(MAPITR_Output$Results[1,2]*100, MAPITR_Output_pValue1, tolerance=1e-8)
	expect_equal(MAPITR_Output$Results[1,3], MAPITR_Output_Est1, tolerance=1e-8)
	expect_equal(MAPITR_Output$Results[1,4], MAPITR_Output_PVE1, tolerance=1e-8)
	expect_equal(MAPITR_Output$Eigenvalues[1,1], MAPITR_Output_Eigen1, tolerance=1e-8)
	expect_equal(MAPITR_Output$Eigenvalues[1,2], MAPITR_Output_Eigen2, tolerance=1e-8)
})

rm(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways, MAPITR_Output, MAPITR_Output_pValue1, MAPITR_Output_Est1, MAPITR_Output_PVE1, MAPITR_Output_Eigen1, MAPITR_Output_Eigen2, envir = .GlobalEnv)
