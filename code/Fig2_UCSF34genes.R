#title: Analysis on UCSF gene signature of 34 genes: GSVA analysis and generating UMAP
#authors: Nayanga Thirimanne

library(clusterProfiler)
library(GSVA)
library(SummarizedExperiment)
library(GSEABase)
library(umap)
library(ggplot2)


#load gene lists
UC <- readxl::read_xlsx("~/UCSF_genepanel_TableS3.xlsx")
UC_E <- UC[UC$Expression_in_meningiomas_with_worse_clinical_outcomes == "Enriched", ] #Upregulated genes
UC_S <- UC[UC$Expression_in_meningiomas_with_worse_clinical_outcomes == "Suppressed", ] #Downregulated genes

#Categorize and label as Enriched and Suppressed genes (= Upregulated and Downregulated genes respectively)
UCg <- UC$Gene %>% as.data.frame()
UCg$term[UCg$. %in% UC_E$Gene] <- "UC_Enriched"
UCg$term[UCg$. %in% UC_S$Gene] <- "UC_Suppressed"

colnames(UCg) <- c("gene", "term")
UCg <- UCg[, c(2,1)] %>% na.omit(.)

UCg_lst=split(UCg, UCg[,1])
UCg_lst=lapply(UCg_lst, function(x) x[,2])

# load batch corrected normalized gene counts for 1298 samples
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")

#calculate GSVA scores
UCg_rec <- GSVA::gsva(log2tpm, UCg_lst, min.sz=10, max.sz=500, verbose = TRUE)

#load UMAP coordinates
umap <- read.csv("~/1298samples_umap2D_041923.csv")

UCg_rec_t <- UCg_rec %>% 
t() %>% 
as.data.frame() %>% 
rownames_to_column(., var = "coordinate_ID")

finaldf_UC <- left_join(umap, UCg_rec_t, by="coordinate_ID")

#plot GSVA scores
ggplot(finaldf_UC,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=UC_Enriched
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("UCSF Enriched genes") +
  scale_colour_gradientn(colours = c(low="navyblue",mid="lightyellow", high="red3"), limits=c(-1,1))
ggsave("UCSF_enrichedgenes_leg.pdf", width = 10, height = 7)
ggsave("UCSF_enrichedgenes_noleg.pdf", width = 7, height = 7)

ggplot(finaldf_UC,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=UC_Suppressed
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("UCSF Suppressed genes") +
  scale_colour_gradientn(colours = c(low="navyblue",mid="lightyellow", high="red3"), limits=c(-1,1)) 
ggsave("UCSF_suppressedgenes_leg.pdf", width = 10, height = 7)
ggsave("UCSF_suppressedgenes_noleg.pdf", width = 7, height = 7)


#ratio of UC_enriched and UC_suppressed
finaldf_UC$UC_Enriched_edit <- finaldf_UC$UC_Enriched + 1
finaldf_UC$UC_Suppressed_edit <- finaldf_UC$UC_Suppressed + 1
finaldf_UC$UC_Ratio <- finaldf_UC$UC_Enriched_edit/finaldf_UC$UC_Suppressed_edit

ggplot(finaldf_UC,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         colour=UC_Ratio
       )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("UCSF Ratio (Enriched score:Suppressed score)") +
  scale_colour_gradientn(colours = c(low="khaki",  high="blue4"))
ggsave("UCSF_Ratio_leg.pdf", width = 10, height = 7)
ggsave("UCSF_Ratio_noleg.pdf", width = 7, height = 7)


# Generating UMAP using UCSF 34 genes and coloring in by metadata
UCgenes <- c(UC_E$Gene, UC_S$Gene)
UC_log <- log2_tpm_adjusted[rownames(log2_tpm_adjusted) %in% UCgenes,]
umap_out <- umap(t(UC_log), random_state = 123, min_dist = 0.1, metric = 'cosine')
umap_2d <- data.frame(umap_out$layout) %>%
  tibble::rownames_to_column("samples")
umap_2d$samples <- gsub("^X", "", umap_2d$samples)
colnames(umap_2d) = c("coordinate_ID","UMAP1_2D", "UMAP2_2D")

sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip = 5)
sample_info <- sample_info[-1, ]
sample_info_2 <- sample_info[match(umap_2d$coordinate_ID, sample_info$coordinate_ID), ] 
finaldf_UC = left_join(umap_2d, sample_info_2,by="coordinate_ID" )

ggplot( finaldf_UC %>% 
          arrange(factor(WHO_GRADE)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(WHO_GRADE),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend+
  ggtitle("WHO Grade") + 
  scale_color_manual(values = c(
    'WHO I' = "orange",
    'WHO II' = "forestgreen",
    'WHO III' = "red",
    "na" = "grey",
    "human YAP1fus" = "blue"
  ))
ggsave("UCSF_34_WHOGRADE_leg.pdf", width = 10, height = 7)
ggsave("UCSF_34_WHOGRADE_noleg.pdf", width = 7, height = 7)

ggplot( finaldf_UC %>% 
          arrange(factor(RECURRENT_TUMOR)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(RECURRENT_TUMOR),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend+
  ggtitle("Recurrent Tumors")  +
  scale_color_manual(values = c(
    "no" = "orange",
    "na" = "grey",
    "yes" = "blue"
  ))
ggsave("UCSF_34_Recurrent_leg.pdf", width = 10, height = 7)
ggsave("UCSF_34_Recurrent_noleg.pdf", width = 7, height = 7)

ggplot( finaldf_UC %>% 
          arrange(factor(CHR22LOSS)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(CHR22LOSS),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend+
  ggtitle("Chr22 Loss")  +
  scale_color_manual(values = c(
    "no" = "orange",
    "na" = "grey",
    "yes" = "blue"
  ))
ggsave("UCSF_34_Chr22loss_leg.pdf", width = 10, height = 7)
ggsave("UCSF_34_Chr22loss_noleg.pdf", width = 7, height = 7)

ggplot( finaldf_UC %>% 
          arrange(factor(UCSF_DNA_METHYLATION_GROUP)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(UCSF_DNA_METHYLATION_GROUP),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend+
  ggtitle("UCSF_DNA_METHYLATION_GROUP")  +
  scale_color_manual(values = c(
    "Hypermitotic" = "red",
    "Immune-enriched" = "forestgreen",
    "Merlin-intact" = "blue",
    "na" = "grey"
  ))
ggsave("UCSF_34_UCSFDNAMethylation_leg.pdf", width = 10, height = 7)
ggsave("UCSF_34_UCSFDNAMethylation_noleg.pdf", width = 7, height = 7)
