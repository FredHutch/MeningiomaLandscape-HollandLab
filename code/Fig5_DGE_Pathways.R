#title: Differential Gene Expression analysis to determine upregulated genes in each cluster
#authors: Nayanga Thirimanne


library(edgeR)

#load cluster IDs from dbscan analysis
clusterID <- readRDS("~/dbscan_1298_cl9_uncl48.rds")
cluster_1 <- clusterID[clusterID$cluster == "1",]
cluster_2 <- clusterID[clusterID$cluster == "2",]
cluster_3 <- clusterID[clusterID$cluster == "3",]
cluster_4 <- clusterID[clusterID$cluster == "4",]
cluster_5 <- clusterID[clusterID$cluster == "5",]
cluster_6 <- clusterID[clusterID$cluster == "6",]
cluster_7 <- clusterID[clusterID$cluster == "7",]
cluster_8 <- clusterID[clusterID$cluster == "8",]
cluster_9 <- clusterID[clusterID$cluster == "9",]

#load batch corrected raw counts
rawcounts1 <- readRDS("/1298_combatseq_rawcounts.rds")

#remove unclustered samples
unclust <- clusterID[clusterID$cluster == "0",]
clusterID <- subset(clusterID, !(clusterID$coordinate_ID %in% unclust$coordinate_ID))
rawcounts1 <- rawcounts1[, !(colnames(rawcounts1) %in% unclust$coordinate_ID)]

#DEG between Cluster A and rest of the samples in the UMAP
#generate sample dataframe
sample <- colnames(rawcounts1)
Adf <- data.frame(sample)
Adf$group[Adf$sample %in% cluster_1$coordinate_ID] <- 'A'
Adf$group[!(Adf$sample %in% cluster_1$coordinate_ID)] <- 'All'
#EdgeR DE analysis: GLM approach (quasi-likelihood F-tests)
#set group
group = factor(Adf$group, levels = c("A","All"))
cols_keep <- Adf$sample
rawcounts <- rawcounts1[, colnames(rawcounts1) %in% cols_keep]
y <- DGEList(counts = rawcounts, genes = rownames(rawcounts), samples = colnames(rawcounts), group = group)
y <- calcNormFactors(y)
#set design
design <- model.matrix(~0+group)
rownames(design) <- colnames(y)
#estimate the NB disperson for the dataset
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
y$tagwise.dispersion
plotBCV(y) 
fit <- glmQLFit(y, design)
my.contrasts1 <- makeContrasts(AvsAll=groupA-groupAll, levels=design)
qlf <- glmQLFTest(fit, contrast = my.contrasts1) 
topTags(qlf)
W <- as.data.frame(topTags(qlf), n=nrow(y))
out <- topTags(qlf, n=Inf)
AvsAll <- as.data.frame(out)
write.csv(AvsAll, "~/Analysis27_DE_1298/DE_dbs_clusterAvsAll_1298_unfilt.csv")

# Above DEG analysis was performed on all the other clusters (e.g. Cluster B vs "All other samples", Cluster C vs "All other samples" etc.)

####
#Plotting pathway enrichment analysis
####

DE <- AvsAll
Upregulated_genes <- DE[DE$logFC <= 0.6 & DE$FDR <= 0.05, ]

#Upregulated genes in each cluster were run on Enrichr.
#Gene Ontology Biological Processes (GOBP) terms enriched in each cluster were downloaded from Enrichr
#Top 15 GOBP terms were plotted as below

plot_title_size = 12
legend_pt_size = 6
axis_text_size = 100
axis_title_size = 15
legend_text_size = 10
spacing = 0.5
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_grid_nolegend <- theme_bw()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank())

theme_grid_legend <- theme_bw()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_text())

En_A <- read.csv("~/Analysis27_DE_1298/Enrichr_output/Enrichr_GO_BP_ClA_edited.csv")
En_A_2 <- top_n(En_A, n = 15, -p.adjust)
En_A_2$Term <- str_replace(En_A_2$Term, " \\s*\\([^\\)]+\\)", "")
En_A_2$Term <- factor(En_A_2$Term, levels = En_A_2$Term)
ggplot(En_A_2, aes(x=-log10(p.adjust), y= reorder(Term, desc(Term)), size=GeneRatio, color=-log10(p.adjust)))+
  geom_point(stat = "identity") +
  scale_color_gradient(low="blue",high="red") +
  ylab("") +
  theme_grid_legend +
  xlab("-log10(p.adjust)")+
  ggtitle("Cluster A")
ggsave("ClusterA_leg.pdf", width = 8, height = 3.5)


### Expression of genes of interest ###

#ggplot2 specifications
plot_title_size = 20
legend_pt_size = 4
axis_text_size = 25
axis_title_size = 25
legend_text_size = 15
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_legend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank() )

#load UMAP 2D coordinates
finaldf <- read.csv("~/1298samples_umap2D_041923.csv", check.names = F)
#load batch corrected normalized gene counts
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")
log2tpm.t <- t(log2tpm)
log2tpm.t <- cbind(rownames(log2tpm.t), data.frame(log2tpm.t, row.names = NULL))
colnames(log2tpm.t)[1] <- "coordinate_ID"


log2tpmt2 <- log2tpm.t[ , c("coordinate_ID", "NF2")]
finaldf_genes <- left_join(finaldf, log2tpmt2, by="coordinate_ID")
my.pallette <- (brewer.pal(9, "Reds"))
ggplot(finaldf_genes,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=NF2
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("NF2 Gene Expression") +
  scale_colour_gradientn(colours = my.pallette)
ggsave("NF2 expression_leg.pdf", width = 10, height = 7)


myPalette <- colorRampPalette((brewer.pal(11, "RdYlBu")))
log2tpmt2 <- log2tpm.t[ , c("coordinate_ID", "PTEN")]
finaldf_genes <- left_join(finaldf_1, log2tpmt2, by="coordinate_ID")
ggplot(finaldf_genes,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=PTEN
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("PTEN") +
  scale_colour_gradientn(colours = myPalette(11))
ggsave("GENEEXP_PTEN_leg.pdf", width = 7, height = 7)

# Similar plots were generated for other genes: HOXD13, HAND2, ROBO1, HOXD12

