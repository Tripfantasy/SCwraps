### Spatial Cell coordinate rotations/flips {Aesthetic}
library(adespatial)

rotate_coords <- function(matrix, angle){
  rotated_matrix <- dummy_matrix
  for(i in 1:nrow(matrix)){
    rotated_matrix[i,] <- rotation(rotated_matrix[i,],angle*pi/180)
  }
  return(rotated_matrix)
}


flip_x <- function(matrix){
  flipped_x <- -matrix$x
  return(flipped_x)
}

flip_y <- function(matrix){
  flipped_y <- -matrix$y
  return(flipped_y)
}


### Basic Seurat Processing {Step 1}
process_objects <- function(wd,pca_dims){
  setwd(wd)
  files <- list.files(pattern = c(".Rds", ".rds"))
  seurat_list <- lapply(files, readRDS)
  for (i in 1:length(seurat_list)) {
    print(paste("Filtering:", i))
    seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^mt-")
    seurat_list[[i]] <- subset(seurat_list[[i]],  percent.mt < 5)
    seurat_list[[i]] <- subset(seurat_list[[i]], subset = nFeature_Spatial < 2500 & nFeature_Spatial > 200)
    print(paste("Preprocessing:", i))
    seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
    seurat_list[[i]] <- ScaleData(seurat_list[[i]], verbose = F)
    seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000)
    seurat_list[[i]] <- SCTransform(seurat_list[[i]], assay = 'Spatial',verbose = FALSE)
    seurat_list[[i]] <- RunPCA(seurat_list[[i]], assay = 'SCT', verbose = FALSE)
    print(paste("Clustering:", i))
    seurat_list[[i]] <- FindNeighbors(seurat_list[[i]], reduction = "pca", dims = 1:pca_dims, verbose = FALSE)
    seurat_list[[i]] <- FindClusters(seurat_list[[i]], verbose = FALSE)
    seurat_list[[i]] <- RunUMAP(seurat_list[[i]], reduction = 'pca', dims = 1:pca_dims, verbose = FALSE)
  }
}

### Split objects {Ctrl v. Exp paradigm}


