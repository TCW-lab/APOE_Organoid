library(Seurat)
library(rhdf5)

setwd('/projectnb/tcwlab/LabMember/akg/APOE3ChPaper/rawData/APOECh_aw')

# read in h5 files using Read10X_h5
a1 <- Read10X_h5('./GSM7729201_alpha_1_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
a2 <- Read10X_h5('./GSM7729202_alpha_2_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
a3 <- Read10X_h5('./GSM7729203_alpha_3_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
a4 <- Read10X_h5('./GSM7729204_alpha_4_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
o1 <- Read10X_h5('./GSM7729205_omega_1_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
o2 <- Read10X_h5('./GSM7729206_omega_2_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
o3 <- Read10X_h5('./GSM7729207_omega_3_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)
o4 <- Read10X_h5('./GSM7729208_omega_4_filtered_feature_bc_matrix.h5', 
                 use.names = TRUE, unique.features = TRUE)

##Create Alpha1-4 & Omega1-4 object and assign it its appropriate genotype values
alpha1 <- CreateSeuratObject(a1)
alpha1@meta.data$APOE <- 'APOE3'
alpha1@meta.data$PSEN <- 'PSEN1'
alpha1@meta.data$E280 <- 'E280A'

alpha2 <- CreateSeuratObject(a2)
alpha2@meta.data$APOE <- 'APOE3Ch'
alpha2@meta.data$PSEN <- 'PSEN1'
alpha2@meta.data$E280 <- 'E280A'

alpha3 <- CreateSeuratObject(a3)
alpha3@meta.data$APOE <- 'APOE3Ch'
alpha3@meta.data$PSEN <- 'PSEN1WT'
alpha3@meta.data$E280 <- '-'

alpha4 <- CreateSeuratObject(a4)
alpha4@meta.data$APOE <- 'APOE3Ch'
alpha4@meta.data$PSEN <- 'PSEN1'
alpha4@meta.data$E280 <- 'E280A'

omega1 <- CreateSeuratObject(o1)
omega1@meta.data$APOE <- 'APOE3'
omega1@meta.data$PSEN <- 'PSEN1WT'
omega1@meta.data$E280 <- '-'

omega1 <- CreateSeuratObject(o2)
omega1@meta.data$APOE <- 'APOE3'
omega1@meta.data$PSEN <- 'PSEN1'
omega1@meta.data$E280 <- 'E280A'

omega3 <- CreateSeuratObject(o3)
omega3@meta.data$APOE <- 'APOE3'
omega3@meta.data$PSEN <- 'PSEN1'
omega3@meta.data$E280 <- 'E280A'

omega4 <- CreateSeuratObject(o4)
omega4@meta.data$APOE <- 'APOE3Ch'
omega4@meta.data$PSEN <- 'PSEN1'
omega4@meta.data$E280 <- 'E280A'

########### merge Seurat Objects ###############################################


organoid0 <- merge(alpha1, y = c(alpha2, alpha3, alpha4), 
                   add.cell.ids = c("alpha1", "alpha2", "alpha3", "alpha4"),
                   project = ""
)





