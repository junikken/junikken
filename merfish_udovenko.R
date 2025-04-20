library(Voyager)
library(SFEData)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater) 
library(ggplot2)
library(patchwork)
library(stringr)
library(spdep)
library(BiocParallel)
library(BiocSingular)
library(gstat)
library(BiocNeighbors)
library(sf)
library(automap)
library(SpatialFeatureExperiment)
options(future.globals.maxSize = 4 * 1024^3)
set.seed(12345)

#Setting dir_use and creating SFE-object
dir_use <- "D:/R-course/Project/uterine_cancer2/"
list.files(dir_use)
mer_sfe <- readVizgen(dir_use, z = 3L, flip = "geometry")

#Getting more info from created SFE object
mer_sfe
#Output:
#class: SpatialFeatureExperiment 
#dim: 550 741749 
#metadata(0):
 # assays(1): counts
#rownames(550): PDK4 CCL26 ... Blank-49 Blank-50
#rowData names(0):
 # colnames(741749): 26985800002100002 26985800002100031 ... 26985800129101535 26985800129101536
#colData names(23): fov volume ... Cellbound1_high_pass sample_id
#reducedDimNames(0):
#  mainExpName: NULL
#altExpNames(0):
#  spatialCoords names(2) : center_x center_y
#imgData names(0):
 #  unit: micron
#Geometries:
 # colGeometries: centroids (POINT), cellSeg (POLYGON) 
#Graphs:
 # sample01: 
summary(mer_sfe)
assayNames(mer_sfe)
head(assay(mer_sfe))

#Visualization
theme_set(theme_bw())
plotGeometry(mer_sfe, colGeometryName = "centroids")
plotGeometry(mer_sfe, colGeometryName = "cellSeg")
plotCellBin2D(mer_sfe, bins = 500, hex = TRUE)


#QC
names(colData(mer_sfe))
system.time(
  print(plotSpatialFeature(mer_sfe, "transcript_count", colGeometryName = "cellSeg",
                      ))
)
#output: user 36.00, system 8.16, elapsed 51.47

system.time({
  print(plotSpatialFeature(mer_sfe, "transcript_count", colGeometryName = "centroids",
                           scattermore = TRUE))
})
#Output: user 2.12, system 0.17, elapsed 2.62

bbox_use <- c(xmin = 3000, xmax = 3500, ymin = 2500, ymax = 3000)
plotSpatialFeature(mer_sfe, "transcript_count", colGeometryName = "cellSeg", bbox = bbox_use)
plotSpatialFeature(mer_sfe, "transcript_count", colGeometryName = "centroids", bbox = bbox_use)
rownames(mer_sfe)
n_genes <- 500

#Plotting distribution of transcript_count divided by number of real genes
colData(mer_sfe)$transcript_count_normed <- mer_sfe$transcript_count / n_genes
plotColDataHistogram(mer_sfe, "transcript_count_normed")

#Plotting distribution of cell volume
plotSpatialFeature(mer_sfe, "volume", colGeometryName = "centroids", scattermore = TRUE)
plotSpatialFeature(mer_sfe, "anisotropy", colGeometryName = "centroids", scattermore = TRUE)

#Relation of transcript detected to cell volume
#plotColData(mer_sfe, x="transcript_count", y="volume", bins=100)
plotColData(mer_sfe, x="volume", y="transcript_count", bins=300)

#Relation volume to anisotropy
plotColData(mer_sfe, x="volume", y="anisotropy", bins=300)


#Transcripts from blank probes
blank_probes <- str_detect(rownames(mer_sfe), "Blank-")
mer_sfe <- addPerCellQCMetrics(mer_sfe, subset = list(blank = blank_probes))
names(colData(mer_sfe))
#Plotting transcript counts from blank probes
plotSpatialFeature(mer_sfe, "subsets_blank_sum", colGeometryName = "centroids", scattermore = TRUE)

#Plotting number of blank probes per cell
plotSpatialFeature(mer_sfe, "subsets_blank_detected", colGeometryName = "centroids", scattermore = TRUE)

#Plotting percentage of blank features per cell
plotSpatialFeature(mer_sfe, "subsets_blank_percent", colGeometryName = "centroids", scattermore = TRUE)

#Plotting blank probes
plotColDataHistogram(mer_sfe, paste0("subsets_blank_", c("sum", "detected", "percent"))) #removal of 2514 rows
mean(mer_sfe$subsets_blank_sum > 0)
#Output: 0,767122 - it means that most cells have at least 1 blank_count

plotColDataHistogram(mer_sfe, "subsets_blank_percent") + scale_x_log10() + annotation_logticks()
#Output: removed 172737 rows containing non-finite outside the scale range (`stat_bin()`).

#Relation of blank percentage to total counts
plotColData(mer_sfe, x="transcript_count", y="subsets_blank_percent", bins=300)
#Output: some cells with high amount of transcript_count have relatively high percent of blank probes 

#Finding outliers
neg_ctrl_outliers <- function(col, mer_sfe, nmads=3, log=FALSE) {
  inds <- colData(mer_sfe)$transcript_count > 0 & colData(mer_sfe)[[col]] > 0
  df <- colData(mer_sfe)[inds, ]
  outliers_inds <- isOutlier(df[[col]], type="higher", nmads = nmads, log=log)
  outliers <- rownames(df)[outliers_inds]
  col2 <- str_remove(col, "^subsets_")
  col2 <- str_remove(col2, "_percent$")
  new_colname <- paste("is", col2, "outlier", sep="_")
  colData(mer_sfe)[[new_colname]] <- colnames(mer_sfe) %in% outliers
  mer_sfe
}
mer_sfe <- neg_ctrl_outliers("subsets_blank_percent", mer_sfe, log=TRUE)
mean(mer_sfe$is_blank_outlier)
#Output: 0.001625887
min(mer_sfe$subsets_blank_percent[mer_sfe$is_blank_outlier])
#Output: 3.508772 (cutoff for outlier)

#Removal of outliers and empty cells
(mer_sfe <- mer_sfe[, !mer_sfe$is_blank_outlier & mer_sfe$transcript_count >0])
#output:
#class: SpatialFeatureExperiment 
#dim: 550 738019 
#metadata(0):
#  assays(1): counts
#rownames(550): PDK4 CCL26 ... Blank-49 Blank-50
#rowData names(1): subsets_blank
#colnames(738019): 26985800003100090 26985800028100439 ... 26985800129101535 26985800129101536
#colData names(31): fov volume ... total is_blank_outlier
#reducedDimNames(0):
#  mainExpName: NULL
#altExpNames(0):
#  spatialCoords names(2) : center_x center_y
#imgData names(1): sample_id
#unit: micron
#Geometries:
#  colGeometries: centroids (POINT), cellSeg (POLYGON) 
#Graphs:
#  sample01: 
#No of cells changed: from 741749 to 738019


#Genes
rowData(mer_sfe)$means <- rowMeans(counts(mer_sfe))
rowData(mer_sfe)$vars <- rowVars(counts(mer_sfe))
rowData(mer_sfe)$blank_probes <- blank_probes
plotRowData(mer_sfe, x="means", y="blank_probes") + scale_y_log10() + annotation_logticks(sides = "b") 
#Overlapping of distrib of some genes with blank

plotRowData(mer_sfe, x="means", y="vars", colour_by = "blank_probes") + geom_abline(slope = 1, intercept = 0, color="green")
  + scale_x_log10() + scale_y_log10() + annotation_logticks() + coord_equal() 

as.data.frame(rowData(mer_sfe)[blank_probes,]) |> ggplot(aes(means, vars)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color="green") +
  scale_x_log10() + scale_y_log10() + annotation_logticks() + coord_equal() 

#Spatial autocorreclation of QC metrics
plotGeometry(mer_sfe, colGeometryName = "cellSeg", bbox = bbox_use)
system.time(colGraph(mer_sfe, "knn5") <- findSpatialNeighbors(mer_sfe, method = "knearneigh",
                                                              dist_type = "idw", k=5, style="W")
            )
#Output: user 72.57, system 0.17, elapsed 73.54

#Computing Moran's I for QC metrics
mer_sfe <- colDataMoransI(mer_sfe, c("transcript_count", "volume"), colGraphName = "knn5")
colFeatureData(mer_sfe)[c("transcript_count", "volume"),]
#Output: moran_sample01 for transcript - 0.243557, volume - 0.163674. k_sample for transcipts - 5.18597, for volume - 6.57556

mer_sfe <- colDataUnivariate(mer_sfe, type="localmoran",
                             features = c("transcript_count", "volume"), colGraphName = "knn5")

plotLocalResult(mer_sfe, "localmoran", c("transcript_count", "volume"), colGeometryName = "centroids",
                scattermore = TRUE, ncol = 2, divergent = TRUE, diverge_center = 0)

#Moran's I
mer_sfe <- logNormCounts(mer_sfe)
system.time(mer_sfe <- runMoransI(mer_sfe, BPPARAM = MulticoreParam(2)))
#Output: user  system elapsed 
        #694.34   30.10  737.53 
plotRowData(mer_sfe, x="moran_sample01", y="blank_probes") + geom_hline(yintercept = 0, linetype = 2)
#Blanks are clustered around 0, most of the genes have positive spatial autocorrelation

#Plotting genes with positive spatial correlation (genes with highest Moran's I)
top_moran <- rownames(mer_sfe)[order(rowData(mer_sfe)$moran_sample01, decreasing = TRUE)[1:16]]
plotSpatialFeature(mer_sfe, top_moran, colGeometryName = "centroids", scattermore = TRUE, ncol = 5)

#Negative Morans
bottom_moran <- rownames(mer_sfe)[order(rowData(mer_sfe)$moran_sample01)[1]]
bottom_abs_moran <- rownames(mer_sfe)[order(abs(rowData(mer_sfe)$moran_sample01))[1]]
plotSpatialFeature(mer_sfe, c(bottom_moran, bottom_abs_moran), colGeometryName = "cellSeg", bbox = bbox_use)


#Spatial autocorrelation (correlogram)
(bins <- st_make_grid(colGeometry(mer_sfe, "centroids"), n = 90, square = FALSE))
#Output: Geometry set for 6862 features 
#Geometry type: POLYGON
#Dimension:     XY
#Bounding box:  xmin: -197.793 ymin: -20.96797 xmax: 12551.22 ymax: 9019.822
#CRS:           NA
#First 5 geometries:
#  POLYGON ((-128.5049 99.04252, -197.793 139.046,...
 #           POLYGON ((-128.5049 339.0635, -197.793 379.067,...
  #                    POLYGON ((-128.5049 579.0845, -197.793 619.088,...
   #                             POLYGON ((-128.5049 819.1055, -197.793 859.109,...
    #                                      POLYGON ((-128.5049 1059.126, -197.793 1099.13,...
df <- cbind(colGeometry(mer_sfe, "centroids"), colData(mer_sfe)[,c("transcript_count", "volume")])
df_binned <- aggregate(df, bins, FUN = mean)
#Removing bins without cells
df_binned <- df_binned[!is.na(df_binned$transcript_count),]
plts <- lapply(c("transcript_count", "volume"), function(f) {
  ggplot(df_binned[,f]) + geom_sf(aes(fill = .data[[f]]), linewidth = 0) + 
    scale_fill_distiller(palette = "Blues", direction = 1) + theme_void()
})
wrap_plots(plts, nrow = 2)

#Computing Moran's I over binned data
nb <- poly2nb(df_binned)
listw <- nb2listw(nb, zero.policy = TRUE)
calculateMoransI(t(as.matrix(st_drop_geometry(df_binned[,c("transcript_count", "volume")]))),
                   listw = listw, zero.policy = TRUE)
#Output: DataFrame with 2 rows and 2 columns
#                     moran         K
#                 <numeric> <numeric>
 # transcript_count  0.606979   14.1201
#   volume          0.400577   62.9936
#Great variance in cells, morans became positive, with increadsing of bin number K value for cell volume also increases and vice versa: with bins=90 K value became 62

#Semivariogram
v <- autofitVariogram(transcript_count ~ 1, as_Spatial(df_binned))
plot(v)

v2 <- variogram(transcript_count ~ 1, data = df_binned, width = 300, cutoff = 4500, map = TRUE)
plot(v2)

#variogram ar specific angles
v3 <- variogram(transcript_count ~ 1, df_binned, alpha = c(30, 90, 150))
v3_model <- fit.variogram(v3, vgm("Ste"))
plot(v3, v3_model)


#PCA
system.time(mer_sfe <- runPCA(mer_sfe, subset_row = !blank_probes, exprs_values = "logcounts", scale = TRUE, BSPARAM = IrlbaParam())
            )
#user  system elapsed 
#236.29    4.88  245.98 

ElbowPlot(mer_sfe)

#Plotting top gene loadings in each PC
plotDimLoadings(mer_sfe)

#Plotting first 4 PCA in space
spatialReducedDim(mer_sfe, "PCA", 4, colGeometryName = "centroids", scattermore = TRUE, divergent = TRUE, diverge_center = 0)

#Zooming into bounding box (components 1, 3, 4)
spatialReducedDim(mer_sfe, "PCA", ncomponents = c(1,3,4), colGeometryName = "cellSeg", bbox = bbox_use, divergent = TRUE, diverge_center = 0)
