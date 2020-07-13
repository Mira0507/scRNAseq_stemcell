#################################### DE analysis ####################################
colData(sce) <- cbind(colData(sce), 
                      clusters) 

library(MAST)

cmatrix <- counts(sce)[!duplicated(rownames(counts(sce))), ]

cdata <- data.frame(wellKey = colnames(counts(sce)),
                    age = colData(sce)$age,
                    cell_type = colData(sce)$lineage,
                    cluster = clusters)

fdata <- data.frame(primerid = rownames(cmatrix))

# Creating SCA (SingleCellAssay) object
sca <- FromMatrix(exprsArray = cmatrix,
                  cData = cdata,
                  fData = fdata,
                  check_sanity = FALSE)


# DEG model by age and cell type
zlm <- zlm(~ age + cell_type + clusters, sca)

# Control Groups: ageold, cell_typeLTHSC, cluster0 
ModelSummary <- summary(zlm)

# fitting the model with contrast 
Fit_Age <- summary(zlm, doLRT = "ageyoung")$datatable
Fit_LTvsST <- summary(zlm, doLRT = "cell_typeSTHSC")$datatable
Fit_LTvsMPP <- summary(zlm, doLRT = "cell_typeMPP")$datatable
Fit_0vs1 <- summary(zlm, doLRT = "clusters1")$datatable


