# Merging and Integration
We started some introductory data entry and analyses for a single visium spatial transcriptomics sample, but what if we want to analyze multiple tissue samples at the same time? This is where the process of merging and integration comes on. The data entry and quality control/preprocessing is very similar to the previous instructions for a single sample, but there are extra steps needed to account for technical differences between the datasets (known as batch effects)!

## Table of Contents

1: Data entry and merging
* 1.1 Data Entry
* 1.2 Merging

2: Quality Control and Preprocessing

3: Cell Clustering and integration
* 3.1: Dimension Reduction and Clustering
* 3.2: Oh No! 
* 3.3: Integration

# Step 1: Data Entry and merging
## 1.1: Data Entry
First step: load all necessary libraries!
```
#Same as before
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

#Some new and improved libraries!
library(glmGamPoi) #this is a newer package that allows for faster SCTransform
library(harmony) #this is the package that will allow us to perform integration and batch corrections
library(UCell) #this is an additional tool for calculating gene signature scores; more on this later
```

The data entry process is the exact same process as before, except we repeat the process for each sample we want to merge:

```
sample <- '"Path/To/Sample/One" #ex: /Users/alfred/Desktop/GSM8633891
object1 <- Load10X_Spatial(data.dir = sample, filename="filtered_feature_bc_matrix.h5")
object1$dataset <- "SampleName1" #assigning unique project name to each object

sample <- '"Path/To/Sample/Two" 
object2 <- Load10X_Spatial(data.dir = sample, filename="filtered_feature_bc_matrix.h5")
object2$dataset <- "SampleName2" #assigning unique project name to each object

#Repeat for each sample that you want to merge

```

The third line (object1$dataset) creates a new metadata column that tracks the dataset of each spatial transcriptomics spot so that it can be traced back to the original dataset during downstream analyses. 

## 1.2: Merging
Next, we want to merge our multiple Seurat objects into a single Seurat object! We can do this using the built in merge() function in Seurat:
```
merged_object <-merge(x = object1, y= list(object2, object3, ...),
                 add.cell.ids = c("SampleName1", "SampleName2", "SampleName3", ...),
                 project = "MergedProject")

#I'm actually not sure if add.cell.ids or project are 100% necessary, I will double check later
```
The result should be a single object (named whatever name you put for "merged_object") containing the merged data of each sample!


# Step 2: Quality Control and Preprocessing (normalization)
Hooray! Now we have a single Seurat object containing the merged data of each sample of interest. Now we proceed with quality control and normalization as we did before (See 1-Introductory-Workflow)

```
#Quick check to see number of cells from each sample!
table(merged_object@meta.data$dataset)

#Calculating mitochondrial and ribosomal DNA percentages
merged_object[["percent.mt"]] <- PercentageFeatureSet(object = merged_object, pattern = "^MT-")
merged_object[["percent.ribo"]] <- PercentageFeatureSet(merged_object, pattern = "^RP[SL]")

#Subsetting based on exclusion criteria
merged_object_subset <- subset(
  merged_object, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 &
  nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)

#Normalization of data using SCTransform
merged_object_subset <- SCTransform(merged_object_subset, assay = "Spatial")

#Quick check to see if the number of viable cells decreased after exclusion
table(merged_object_subset@meta.data$dataset)

```

# Step 3: Cell Clustering and integration
As you may remember, the next step after quality control and data preprocessing is to run clustering algorithms on our data to identify populations of cells with similar gene expression data. This allows us to better isolate specific cell types for downstream analyses.

## 3.1: Dimension Reduction and clustering
The instructions for this step are similar to before as well:
```
#Standard workflow for running PCA clustering and UMAP
merged_object_subset <- RunPCA(merged_object_subset, assay = "SCT", reduction.name = "pca.SCT")
merged_object_subset <- FindNeighbors(merged_object_subset, assay = "SCT", reduction = "pca.SCT", dims = 1:15)
merged_object_subset <- FindClusters(merged_object_subset, cluster.name = "seurat_cluster.SCT", resolution = 0.5)
merged_object_subset <- RunUMAP(merged_object_subset, reduction = "pca.SCT", reduction.name = "umap.SCT", return.model = T, dims = 1:15)

#Display UMAP plot, using DimPLot()
DimPlot(SeuratObject_subset, reduction = "umap.SCT", label = TRUE)  #using the reduction from the last UMAP line
```

You should see the data clustered into distinct cell populations. Looks great! ... right?
![Alt text](https://raw.githubusercontent.com/alfalfacow/Spatial-Transcriptomics/main/Images/fake-umap.png)

## 3.2: Oh No!
As you can tell by the title of this section, there is a big problem!!! As mentioned in the introduction of this page, merging multiple datasets introduces batch effects, in which different samples are "biased" in a sense due to the unique technical environments that went into collecting each of these samples and the data involved.

To see the effect of batch effects, we can display the UMAP plot, with the spots labeled by dataset instead of by the UMAP clusters:
```
#Use the Idents() function to change the main "identity" label of the cells, and thus the labels on the UMAP
Idents(merged_object_subset) <- "dataset" #sets default labeling on plots to be datasets, USE THIS FOR THIS STEP
Idents(merged_object_subset) <- "seurat_cluster.SCT" #sets default labeling on plots to be the clusters, USE THIS TO REVERT BACK TO CLUSTERS

#Display UMAP plot again, after changing labelling of cells by dataset (which we had defined earlier to be the sample name)
DimPlot(PNI, reduction = "umap.SCT", label = TRUE)
```

The resulting UMAP plot should have a major problem: the cells seem to be clustered together based on dataset type! Here is an example:
![Alt text](https://raw.githubusercontent.com/alfalfacow/Spatial-Transcriptomics/main/Images/bad-umap.png)
