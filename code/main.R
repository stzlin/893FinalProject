library('R.matlab')
setwd("/Users/stephanie/UNC-chapel hill/Spring2021/STOR893/FinalProject/data/")
filename = system.file('SC/HCP_cortical_DesikanAtlas_SC.mat', package = 'R.matlab')

path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path("SC", "HCP_cortical_DesikanAtlas_SC.mat")
SC <- readMat(pathname)

pathname <- file.path("FC", "HCP_cortical_DesikanAtlas_FC.mat")
FC <- readMat(pathname)

pathname1 <- file.path("TNPCA_Result", "TNPCA_Coeff_HCP_Structural_Connectome.mat")
pathname2 <- file.path("TNPCA_Result", "TNPCA_Coeff_HCP_Functional_Connectome.mat")
TNPCA_Structural <- readMat(pathname1)
TNPCA_Functional <- readMat(pathname2)

pathname <- file.path("traits/175traits", "HCP_175Traits.mat")
traits<-readMat(pathname)
