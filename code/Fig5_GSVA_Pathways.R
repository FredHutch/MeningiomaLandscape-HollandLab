
#title: GSVA analysis on biological pathways 
#authors: Nayanga Thirimanne

library(ggplot2)
library(clusterProfiler)
library(GSVA)
library(SummarizedExperiment)
library(GSEABase)
library(edgeR)

setwd("/Meningioma_bulkrnaseq_datasets/msigdb_v7.2_GMTs/")
biocarta_db <- read.gmt("c2.cp.biocarta.v7.2.symbols.gmt")
gobp_db <- read.gmt("c5.go.bp.v7.2.symbols.gmt")


biocarta_lst=split(biocarta_db, biocarta_db[,1])
biocarta_lst=lapply(biocarta_lst, function(x) x[,2])


gobp_lst=split(gobp_db, gobp_db[,1])
gobp_lst=lapply(gobp_lst, function(x) x[,2])

#load normalized gene counts
my_mat <- readRDS("~/1298samples_combatseq_log2tpm.rds")

biocarta_res <- GSVA::gsva(my_mat, biocarta_lst, min.sz=10, max.sz=500, verbose=TRUE)
gobp_res <- GSVA::gsva(my_mat, gobp_lst, min.sz=10, max.sz=500, verbose=TRUE)

score=rbind(biocarta_res, gobp_res)
saveRDS(score, "~/1298_gsvascored_biocarta.gobp.rds")

## Color in UMAP by GSVA scores for specific pathways ##

#load UMAP coordinates
umap <- read.csv("~/1298samples_umap2D_041923.csv")
#load GSVA scores
gsva_1271 <- readRDS("~/1298_gsvascores_gobp.rds")
gsva_1271.t <- t(gsva_1271)
gsva_1271.t <- cbind(rownames(gsva_1271.t), data.frame(gsva_1271.t, row.names = NULL))
colnames(gsva_1271.t)[1] <- "coordinate_ID"

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#GOBP: Cell Cycle DNA Replication 
gsva_1271_2 <- gsva_1271.t[ , c("coordinate_ID", "GO_CELL_CYCLE_DNA_REPLICATION")]
finaldf_genes <- left_join(finaldf, gsva_1271_2, by="coordinate_ID")
ggplot(finaldf_genes %>% arrange(GO_CELL_CYCLE_DNA_REPLICATION),
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=GO_CELL_CYCLE_DNA_REPLICATION
       )
) + 
  geom_point(size=1.5) +
  theme_nolegend +
  ggtitle("CELL CYCLE DNA REPLICATION") +
  scale_colour_gradientn(colours = myPalette(11), limits=c(-1,1))
ggsave("GO_CELLCYCLE.pdf", width = 7, height = 7)
ggsave("GO_CELLCYCLE_noleg.pdf", width = 7, height = 7)

#GOBP: Neurotransmitter secretion
gsva_1271_2 <- gsva_1271.t[ , c("coordinate_ID", "GO_NEUROTRANSMITTER_SECRETION")]
finaldf_genes <- left_join(finaldf, gsva_1271_2, by="coordinate_ID")
ggplot(finaldf_genes,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=GO_NEUROTRANSMITTER_SECRETION
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("NEUROTRANSMITTER SECRETION") +
  scale_colour_gradientn(colours = myPalette(11), limits=c(-1,1))
ggsave("GO_NEUROTRANSMITTER_SECRETION.pdf", width = 7, height = 7)
ggsave("GO_NEUROTRANSMITTER_SECRETION_noleg.pdf", width = 7, height = 7)

#GOBP:Antigen Processing and Presentation
gsva_1271_2 <- gsva_1271.t[ , c("coordinate_ID", "GO_ANTIGEN_PROCESSING_AND_PRESENTATION")]
finaldf_genes <- left_join(finaldf, gsva_1271_2, by="coordinate_ID")
ggplot(finaldf_genes %>% arrange(GO_ANTIGEN_PROCESSING_AND_PRESENTATION),
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=GO_ANTIGEN_PROCESSING_AND_PRESENTATION
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("ANTIGEN PROCESSING_AND PRESENTATION") +
  scale_colour_gradientn(colours = myPalette(11), limits=c(-1,1))
ggsave("GO_ANTIGEN_PROCESSING_AND_PRESENTATION.pdf", width = 7, height = 7)
ggsave("GO_ANTIGEN_PROCESSING_AND_PRESENTATION_noleg.pdf", width = 7, height = 7)

#GOBP: Vasculogenesis
gsva_1271_2 <- gsva_1271.t[ , c("coordinate_ID", "GO_VASCULOGENESIS")]
finaldf_genes <- left_join(finaldf, gsva_1271_2, by="coordinate_ID")
ggplot(finaldf_genes,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=GO_VASCULOGENESIS
       )
) + 
  geom_point(size=1.5) +
  theme_nolegend +
  ggtitle("VASCULOGENESIS") +
  scale_colour_gradientn(colours = myPalette(11), limits=c(-1,1))
ggsave("GO_VASCULOGENESIS.pdf", width = 7, height = 7)
ggsave("GO_VASCULOGENESIS_noleg.pdf", width = 7, height = 7)

