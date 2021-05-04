library(tidyverse)
library(dplyr)
library(Seurat)
reticulate::py_install(packages = 'umap-learn')


#set working directory
setwd("/projectnb/bf528/users/swiss_cheese2/project_4/analyst")

#sample seurat object data 
cells <- readRDS("GSM2230760_seurat.rda")

#4.1 
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
cells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# Filter on only significant markers:
cells.markers.padj <- cells.markers[cells.markers$p_val_adj<0.05,]
cells.markers.padj$gene

# Save to csv file:
write.csv(cells.markers, file="diff_genes.csv")

## finding cell types found in paper that are present in our data set
# Epsilon, vascular, cytotoxic and mast cells are absent 


#4.2 Label marker genes found in clusters -- marker genes are from paper 
#GCG alpha 2,4,8,10,11 
#INS beta 1, 2, 6
#SST delta 0, 3 
#PPY gamma 0 
#GHRL Epsilon NONE 
#KRT19 Ductal 5,10 
#CPA1 Acinar 3, 11 
#PDGFRA fibroblasts 11 
#PDGFRB fibroblasts 9 
#VWF, PECAM1,CD34 NONE 
#CD68 Macrophage 10, 11, 12 
#TPSAB1, KIT, CPA3 Mast NONE



top5 <-cells.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top30 <- cells.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

cluster0.marker <- top5[top5$cluster==0,]
cluster1.marker <- top5[top5$cluster==1,]
cluster2.marker <- top5[top5$cluster==2,]
cluster3.marker <- top5[top5$cluster==3,]
cluster4.marker <- top5[top5$cluster==4,]
cluster5.marker <- top5[top5$cluster==5,]
cluster6.marker <- top5[top5$cluster==6,]
cluster7.marker <- top30[top30$cluster==7,]
cluster8.marker <- top30[top30$cluster==8,]
cluster9.marker <- top5[top5$cluster==9,]
cluster10.marker <- top5[top5$cluster==10,]
cluster11.marker <- top30[top30$cluster==11,]
cluster12.marker <- top30[top30$cluster==12,]


new.cluster.ids <- c("Delta", "Beta", "Alpha_1", "Acinar", "Alpha_2","Ductal","Beta_2", 
                     "Macrophage", "Alpha_3","Stellate",
                     "Cytotoxic T","Cytotoxic/Acinar", "Macrophage")

names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
Umap <- DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 1,
                label.size = 5) + NoLegend()
Umap

top1 <- cells.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#heatmap 
top5 <- cells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

map <- DoHeatmap(cells, features=top1$gene, size=5, angle=75)

# Save heatmap:
png(filename="heatmap_top5.png", width=1200, height=800)
map
dev.off()




