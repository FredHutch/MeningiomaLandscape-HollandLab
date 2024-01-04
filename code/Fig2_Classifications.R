#title: Color in UMAP by tumor classification methods
#authors: Nayanga Thirimanne


library(ggplot2)

#read count file with umap coordinates
finaldf <- read.csv("~/1298samples_umap2D_041923.csv", check.names = FALSE)


# UCSF DNA Methylation 
ggplot( finaldf %>% 
          arrange(desc(factor(UCSF_DNA_METHYLATION_GROUP))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(UCSF_DNA_METHYLATION_GROUP),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("UCSF DNA METHYLATION GROUP") +
  scale_color_manual(values = c("Hypermitotic" = "red",
                                "Immune-enriched" = "forestgreen",
                                "Merlin-intact" = "blue",
                                "na" = "grey"))

# Baylor RNA Classifier
ggplot( finaldf %>% 
          arrange(desc(factor(BAYLOR_RNA_CLASSIFIER))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(BAYLOR_RNA_CLASSIFIER),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("BAYLOR RNA CLASSIFIER") +
  scale_color_manual(values = c("A/NF2wt_ben" = "blue",
                                "B/NF2loss_int" = "forestgreen",
                                "C/NF2loss_mal" = "red",
                                "na" = "grey"))

# Toronto Methyl Profile Grade
ggplot( finaldf %>% 
          arrange(desc(factor(TORONTO_METHYL_PROFILE_GRADE_COCA))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(TORONTO_METHYL_PROFILE_GRADE_COCA),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("TORONTO METHYL PROFILE GRADE") +
  scale_color_manual(values = c("MG1/Immunogenic" = "chocolate1",
                                "MG2/Benign_NF2wt" = "blue",
                                "MG3/Hypermetabolic" = "green4",
                                "MG4/Proliferative" = "red",
                                "na" = "grey"))

# Ferreira Grade
ggplot( finaldf %>% 
          arrange(desc(factor(FERREIRA_GRADE))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(FERREIRA_GRADE),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("FERREIR GRADE") +
  scale_color_manual(values = c("1" = "gold",
                                "1.5" = "blue",
                                "2" = "forestgreen",
                                "3" = "red",
                                "na" = "grey"))

# Simpson Grade
ggplot( finaldf %>% 
          arrange(desc(factor(SIMPSON_GRADE))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(SIMPSON_GRADE),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("SIMPSON GRADE") +
  scale_color_manual(values = c("1" = "gold",
                                "2" = "blue",
                                "3" = "green3",
                                "4" = "red",
                                "na" = "grey"))

# Sahm Classification
ggplot( finaldf %>% 
          arrange(factor(SAHM_METH_CLASS)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(SAHM_METH_CLASS),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("Sahm classification") +
  scale_color_manual(values = c("human YAP1fus" = "black",
                                "NF2 WT Ben-2" = "red",
                                "NF2 Mut Ben-1" = "blue",
                                "NF2 Mut Mal" = "green",
                                "na" = "grey"))

