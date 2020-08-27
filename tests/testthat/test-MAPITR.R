context("Tests for MAPITR.R") 

data(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways, MAPITR_TestData_PCs)

assign("MAPITR_TestData_Genotypes", MAPITR_TestData_Genotypes, envir = .GlobalEnv)
assign("MAPITR_TestData_Phenotype", MAPITR_TestData_Phenotype, envir = .GlobalEnv)
assign("MAPITR_TestData_Pathways", MAPITR_TestData_Pathways, envir = .GlobalEnv)
assign("MAPITR_TestData_PCs", MAPITR_TestData_PCs, envir = .GlobalEnv)
assign("MAPITR_Output1", MAPITR(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways), envir = .GlobalEnv)
assign("MAPITR_Output1_pValue1", 0.06793231, envir = .GlobalEnv)
assign("MAPITR_Output1_Est1", 1.0684080, envir = .GlobalEnv)
assign("MAPITR_Output1_PVE1", 1.1281354, envir = .GlobalEnv)
assign("MAPITR_Output1_Eigen1", -0.02115170, envir = .GlobalEnv)
assign("MAPITR_Output1_Eigen2", -0.01870068, envir = .GlobalEnv)
assign("MAPITR_Output2", MAPITR(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways, Covariates=MAPITR_TestData_PCs), envir = .GlobalEnv)
assign("MAPITR_Output2_pValue1", 0.01282748, envir = .GlobalEnv)
assign("MAPITR_Output2_Est1", 1.1540813, envir = .GlobalEnv)
assign("MAPITR_Output2_PVE1", 1.2908323, envir = .GlobalEnv)
assign("MAPITR_Output2_Eigen1", -0.01983171, envir = .GlobalEnv)
assign("MAPITR_Output2_Eigen2", -0.01788879, envir = .GlobalEnv)

test_that("MAPITR runs the main MAPITR function", {
	expect_equal(unname(MAPITR_Output1$Results[1,2]*100), MAPITR_Output1_pValue1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output1$Results[1,3]), MAPITR_Output1_Est1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output1$Results[1,4]), MAPITR_Output1_PVE1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output1$Eigenvalues[1,1]), MAPITR_Output1_Eigen1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output1$Eigenvalues[1,2]), MAPITR_Output1_Eigen2, tolerance=1e-4)
})

test_that("MAPITR runs the main MAPITR function with covariate", {
	expect_equal(unname(MAPITR_Output2$Results[1,2]*100), MAPITR_Output2_pValue1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output2$Results[1,3]), MAPITR_Output2_Est1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output2$Results[1,4]), MAPITR_Output2_PVE1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output2$Eigenvalues[1,1]), MAPITR_Output2_Eigen1, tolerance=1e-4)
	expect_equal(unname(MAPITR_Output2$Eigenvalues[1,2]), MAPITR_Output2_Eigen2, tolerance=1e-4)
})

rm(MAPITR_TestData_Genotypes, MAPITR_TestData_Phenotype, MAPITR_TestData_Pathways, MAPITR_TestData_PCs, MAPITR_Output1, MAPITR_Output1_pValue1, MAPITR_Output1_Est1, MAPITR_Output1_PVE1, MAPITR_Output1_Eigen1, MAPITR_Output1_Eigen2, MAPITR_Output2, MAPITR_Output2_pValue1, MAPITR_Output2_Est1, MAPITR_Output2_PVE1, MAPITR_Output2_Eigen1, MAPITR_Output2_Eigen2, envir = .GlobalEnv)
