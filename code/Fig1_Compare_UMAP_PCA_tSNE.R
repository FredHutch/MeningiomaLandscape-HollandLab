#title: Compare dimension reduction techniques (UMAP, PCA and tSNE) with batchcorrected and non-batchcorrected datasets
#authors: Nayanga Thirimanne

library(ggplot2)
library(umap)
library(Rtsne)
library(dplyr)
library(readxl)
library(tibble)
library(factoextra)

# ggplot2 specifications 
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
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank() )

theme_nolegend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank() )

theme_grid_nolegend <- theme_bw()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank())

#################################
#Compare UMAP, PCA and tSNE using normalized non-batchcorrected datasets
#################################

#load count file
log2_tpm_noBC <- readRDS("~/1298_normalized_nonbatchcorrected.rds")

#generate UMAP coordinates
umap_out <- umap(t(log2_tpm_noBC), random_state = 123, min_dist = 0.1, metric = 'cosine')
umap_2d <- data.frame(umap_out$layout) %>%
  tibble::rownames_to_column("samples")
umap_2d$samples <- gsub("^X", "", umap_2d$samples)
colnames(umap_2d) = c("coordinate_ID","UMAP1_2D", "UMAP2_2D")
#combine with metadata file
sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip = 5)
sample_info <- sample_info[-1, ]
sample_info_2 <- sample_info[match(umap_2d$coordinate_ID, sample_info$coordinate_ID), ] 
finalDf2 = left_join(umap_2d, sample_info_2,by="coordinate_ID" )

#color in UMAP by datasets
ggplot( finalDf2 %>% 
          arrange(factor(DATASET)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(DATASET),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "chocolate",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "red4",
                                "UW/FHCC" = "red"))
ggsave("DATASET_woBC_leg.pdf", width = 10, height = 7)
ggsave("DATASET__woBC_noleg.pdf", width = 7, height = 7)

#PCA
pca_matrix <- as.matrix(log2_tpm_noBC) %>% t()
sample_pca <- prcomp(pca_matrix, scale. = FALSE)
pc_scores <- sample_pca$x
variances <- (100*((sample_pca$sdev)^2))/sum((sample_pca$sdev)^2)
variances
pc_scores %>% 
  as_tibble(rownames = "coordinate_ID") %>% 
  full_join(sample_info_2, by="coordinate_ID" ) %>%
  ggplot(aes(x=PC1, y=PC2, colour=factor(DATASET))) +
  geom_point() +
  theme_bw() +
  xlab("PC1 (37%)")+
  ylab("PC2 (8%)") +
  theme_grid_legend +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "chocolate",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "red4",
                                "UW/FHCC" = "red"))

ggsave("PCA_woBC_leg.pdf", width = 7, height = 7)
ggsave("PCA_woBC_noleg.pdf", width = 7, height = 7)

#tSNE
tpm_2 <- t(log2_tpm_noBC)
tsne_result <- Rtsne(tpm_2, perplexity = 30, check_duplicates = FALSE)
tsne <- tsne_result$Y
rownames(tsne) <- rownames(tpm_2)
t2 <- tsne %>% 
  as_tibble(rownames = "coordinate_ID") %>% 
  full_join(sample_info_2, by="coordinate_ID" )
ggplot(t2, aes(x=V1, y=V2, colour=factor(DATASET))) +
  geom_point() +
  theme_bw() +
  xlab("tSNE1") +
  ylab("tSNE2") +
  theme_grid_nolegend +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "chocolate",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "red4",
                                "UW/FHCC" = "red"))
ggsave("tSNE_woBC_leg.pdf", width = 7, height = 7)
ggsave("tSNE_woBC_noleg.pdf", width = 7, height = 7)

#################################
#Compare UMAP, PCA and tSNE using normalized batchcorrected datasets
#################################

#load count file
log2_tpm_adjusted <- readRDS("~/Output_1298_April2023/1298samples_combatseq_log2tpm.rds")

#UMAP
umap_out <- umap(t(log2_tpm_adjusted), random_state = 123, min_dist = 0.1, metric = 'cosine')
umap_2d <- data.frame(umap_out$layout) %>%
  tibble::rownames_to_column("samples")
umap_2d$samples <- gsub("^X", "", umap_2d$samples)
colnames(umap_2d) = c("coordinate_ID","UMAP1_2D", "UMAP2_2D")
#Combine with metadata
sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip = 5)
sample_info <- sample_info[-1, ]
sample_info_2 <- sample_info[match(umap_2d$coordinate_ID, sample_info$coordinate_ID), ] 
finaldf = left_join(umap_2d, sample_info_2,by="coordinate_ID" )

#color in UMAP by datasets
ggplot( finaldf %>% 
          arrange(factor(DATASET)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(DATASET),
        )
) + 
  geom_point(size=1.5) +
  theme_legend+
  ggtitle("Datasets") + 
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "brown",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "brown",
                                "UW/FHCC" = "red"))
ggsave("DATASET_leg.pdf", width = 10, height = 7)

#PCA
pca_matrix <- as.matrix(log2_tpm_adjusted) %>% t()
sample_pca <- prcomp(pca_matrix, scale. = FALSE)
pc_scores <- sample_pca$x
variances <- (100*((sample_pca$sdev)^2))/sum((sample_pca$sdev)^2)
variances
pc_scores %>% 
  as_tibble(rownames = "coordinate_ID") %>% 
  full_join(sample_info_2, by="coordinate_ID" ) %>%
  ggplot(aes(x=PC1, y=PC2, colour=factor(DATASET))) +
  geom_point() +
  theme_bw() +
  xlab("PC1 (27%)")+
  ylab("PC2 (6%)") +
  theme_grid_legend +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "chocolate",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "red4",
                                "UW/FHCC" = "red"))

ggsave("PCA_leg.pdf", width = 7, height = 7)

#tSNE
tpm_2 <- t(log2_tpm_adjusted)
tsne_result <- Rtsne(tpm_2, perplexity = 30, check_duplicates = FALSE)
tsne <- tsne_result$Y
rownames(tsne) <- rownames(tpm_2)
t2 <- tsne %>% 
  as_tibble(rownames = "coordinate_ID") %>% 
  full_join(sample_info_2, by="coordinate_ID" )
ggplot(t2, aes(x=V1, y=V2, colour=factor(DATASET))) +
  geom_point() +
  theme_bw() +
  xlab("tSNE1") +
  ylab("tSNE2") +
  theme_grid_nolegend +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c("Baylor" = "blue",
                                "CAVATICA" = "grey",
                                "CHOP" = "forestgreen",
                                "Yale" = "chocolate",
                                "Heidelberg" = "purple",
                                "HKU/UCSF" = "deepskyblue",
                                "Palacky" = "black",
                                "UCSD" = "violet",
                                "UCSF(2018)" = "cyan",
                                "UCSF(2020)" = "gold", 
                                "UCSF(2022)" = "orange",
                                "UoT" = "red4",
                                "UW/FHCC" = "red"))
ggsave("tSNE_leg.pdf", width = 7, height = 7)

#color in PCA and tSNE by WHO Grade
pc_scores %>% 
  as_tibble(rownames = "coordinate_ID") %>% 
  full_join(sample_info_2, by="coordinate_ID" ) %>%
  ggplot(aes(x=PC1, y=PC2, colour=factor(WHO_GRADE))) +
  geom_point() +
  theme_bw() +
  xlab("PC1 (27%)")+
  ylab("PC2 (6%)") +
  theme_grid_legend + 
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c(
    'WHO 1' = "orange",
    'WHO 2' = "forestgreen",
    'WHO 3' = "red",
    "na" = "grey"
  ))  
ggsave("PCA_WHO_BC_leg.pdf", width = 7, height = 7)

ggplot(t2, aes(x=V1, y=V2, colour=factor(WHO_GRADE))) +
  geom_point() +
  theme_bw() +
  xlab("tSNE1") +
  ylab("tSNE2") +
  theme_grid_legend +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c(
    'WHO 1' = "orange",
    'WHO 2' = "forestgreen",
    'WHO 3' = "red",
    "na" = "grey"
  )) 
ggsave("tSNE_WHO_BC_leg.pdf", width = 7, height = 7)

#################

