#installing and calling packages 
BiocManager::install("tximport")
BiocManager::install("SeqGSEA")
BiocManager::install('fishpond')
BiocManager::install('EnsDb.Hsapiens.v79')
install.packages("Seurat")
install.packages("Matrix")
install.packages("tidyverse")
library('tidyverse')
library('Matrix')
library('Seurat')
library('tximport')
library('SeqGSEA')
library('fishpond')
library('EnsDb.Hsapiens.v79')

x <- file.path("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz")  #defining file path
txi<- tximport(x, type="alevin")  #reading in alevin counts using txi
t<- txi$counts

#creating Seurat Object
panc_cells <- CreateSeuratObject(counts = txi$counts, project = "panc_cells", min.cells = 3, min.features = 200)
panc_cells
head(panc_cells)

#mapping gene symbols
ids <- rownames(t)
View(ids)
ids <- sub("[.][0-9]*", "", ids)
dim(ids)
rownames(t) <- ids
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

t <- t[rownames(t) %in% geneIDs$GENEID,]
rownames(t) <- geneIDs$SYMBOL

panc_cells[["percent.mt"]] <- PercentageFeatureSet(panc_cells, pattern = "^MT-") #calculates percentages of counts and does quality check
head(panc_cells@meta.data, 5)

VlnPlot(panc_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2) #visualises as violin plot

plot1 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filtering low quality cells 
panc_cells <- subset(panc_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize the data
panc_cells <- NormalizeData(panc_cells)

#identify highly variable features 
panc_cells <- FindVariableFeatures(panc_cells, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(panc_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the data
all.genes <- rownames(panc_cells)
panc_cells <- ScaleData(pbmc, features = all.genes)

#perform linear dimensional reduction
panc_cells <- RunPCA(panc_cells, features = VariableFeatures(object = panc_cells))
print(panc_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(panc_cells, dims = 1:2, reduction = "pca")
DimPlot(panc_cells, reduction = "pca")

#testing for association between observed and latent variables
panc_cells <- JackStraw(panc_cells, num.replicate = 100)
panc_cells <- ScoreJackStraw(panc_cells, dims = 1:20)
JackStrawPlot(panc_cells, dims = 1:15)
ElbowPlot(panc_cells)

#clustering the cells
panc_cells <- FindNeighbors(panc_cells, dims = 1:10)
panc_cells <- FindClusters(panc_cells, resolution = 0.5)

#run non-linear reduction
panc_cells <- RunUMAP(panc_cells, dims = 1:10)
DimPlot(panc_cells, reduction = "umap")
