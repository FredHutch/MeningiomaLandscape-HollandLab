#title: Visualizing gene fusion calls made using Arriba
#authors: Nayanga Thirimanne

library(ggplot2)

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


theme_nolegend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank() )

# Calculating fusion burden 
# Only the high confidence fusions, with atleast one coding gene partner that recur at least twice within dataset were taken into account
FusJ <- read.csv("~/at_least_one_fusion_gene_are_coding_1249x2896.csv", check.names = FALSE, row.names = 1)
FusJ2 <- colSums(FusJ) %>% as.data.frame()
FusJ2 <- rownames_to_column(FusJ2, "Fusions")
colnames(FusJ2)[2] <- "Freq"
FusJ3 <- FusJ2[order(FusJ2$Freq, decreasing = TRUE),]
FusJ4 <- FusJ3[FusJ3$Freq >=2, ]

#FusionBurden - calculated
#only consider fusions that are present at least twice
FusJ.1 <- FusJ[, colnames(FusJ) %in% FusJ4$Fusions]
write.csv(FusJ.1, "~/HighConf_ProteinCoding_Rec2_171fusions.csv")
FusJ5 <- rowSums(FusJ.1) %>% as.data.frame()
FusJ5 <- rownames_to_column(FusJ5, "coordinate_ID")
colnames(FusJ5)[2] <- "FusionBurden"
finaldf_1 <- finaldf[, 1:4]
finaldf_f <- left_join(finaldf_1, FusJ5, by="coordinate_ID" )

finaldf_f$FusionBurden <- as.numeric(finaldf_f$FusionBurden)
finaldf_f$FusionBurden_category <- cut(finaldf_f$FusionBurden, breaks = c(-1, 0, 1, 3, 5, 7, 10, 20 ),
                                       labels = c("0", "1", "2-3", "4-5", "6-7", "8-10", "10<"))
finaldf_f$FusionBurden_category <- as.character(finaldf_f$FusionBurden_category)
finaldf_f$FusionBurden_category[is.na(finaldf_f$FusionBurden_category)] <- "na"
finaldf_f$FusionBurden_category <- as.factor(finaldf_f$FusionBurden_category)

One <- subset(finaldf_f, finaldf_f$FusionBurden_category == "1")
Two1 <- subset(finaldf_f, finaldf_f$FusionBurden_category == "2-3")
Two <- subset(finaldf_f, finaldf_f$FusionBurden_category == "4-5")
Three <- subset(finaldf_f, finaldf_f$FusionBurden_category == "6-7")
Four <- subset(finaldf_f, finaldf_f$FusionBurden_category == "8-10")
Five <- subset(finaldf_f, finaldf_f$FusionBurden_category == "10<")
Six <- subset(finaldf_f, finaldf_f$FusionBurden_category == "na")
Seven <- subset(finaldf_f, finaldf_f$FusionBurden_category == "0")

finaldf_f$FusionBurden_category <- factor(finaldf_f$FusionBurden_category, levels= c("0", "1-3", "4-5", "6-7", "8-10", "10<"))

ggplot(finaldf_f,
       aes(
         x = UMAP1_2D, 
         y = UMAP2_2D,
         color = factor(FusionBurden_category),
       )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("Fusion burden") + 
  geom_point(data = Six, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Seven, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Five, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = One, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Two, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Two1, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Four, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  geom_point(data = Three, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(FusionBurden_category))) +
  scale_color_manual(values = c(
    '1' = "cyan3",
    '2-3' = "darkcyan",
    '4-5' = "orchid1",
    '6-7' = "red",
    '8-10' = "blue",
    '10<' = "black",
    "na"  = "grey80", 
    "0" = "grey80"
  ))
ggsave("FusionBurden_leg_new_rec2.pdf", width = 10, height = 7)
ggsave("FusionBurden_noleg_new_rec2.pdf", width = 7, height = 7)

##
#Plotting fusions of interest
##

# TRPM3(105635),RP11-563H8.2(37062)--TRPM3
Fus3 <- rownames_to_column(FusJ, var = "coordinate_ID")
finaldf_g <- left_join(finaldf_f, Fus3, by="coordinate_ID" ) #finaldf_f from above code

finaldf_g$`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`[is.na(finaldf_g$`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`)] <- "na"
FF <- subset(finaldf_g, finaldf_g$`TRPM3(105635),RP11-563H8.2(37062)-TRPM3` == 1) 
FF2 <- subset(finaldf_g, finaldf_g$`TRPM3(105635),RP11-563H8.2(37062)-TRPM3` == 0) 
FF3 <- subset(finaldf_g, finaldf_g$`TRPM3(105635),RP11-563H8.2(37062)-TRPM3` == "na")
ggplot( finaldf_g %>% 
          arrange(factor(`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`),
        )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`))) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`TRPM3(105635),RP11-563H8.2(37062)-TRPM3`))) +
  theme_legend + 
  ggtitle("TRPM3(105635),RP11-563H8.2(37062)-TRPM3") + 
  scale_color_manual(values = c(
    '0' = "grey",
    '1' = "red",
    '2' = "blue",
    "3" = "green",
    "4" = "violet", 
    "na"  = "grey60"
  ))

ggsave("Fusion_TRPM3_leg.pdf", width = 10, height = 7)
ggsave("Fusion_TRPM3_noleg.pdf", width = 7, height = 7)

# PARD6B(19948),BCAS4(18151)--BCAS4
finaldf_g$`PARD6B(19948),BCAS4(18151)-BCAS4`[is.na(finaldf_g$`PARD6B(19948),BCAS4(18151)-BCAS4`)] <- "na"
FF <- subset(finaldf_g, finaldf_g$`PARD6B(19948),BCAS4(18151)-BCAS4` == 1) 
FF2 <- subset(finaldf_g, finaldf_g$`PARD6B(19948),BCAS4(18151)-BCAS4` == 0) 
FF3 <- subset(finaldf_g, finaldf_g$`PARD6B(19948),BCAS4(18151)-BCAS4` == "na")
ggplot( finaldf_g %>% 
          arrange(factor(`PARD6B(19948),BCAS4(18151)-BCAS4`)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(`PARD6B(19948),BCAS4(18151)-BCAS4`),
        )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`PARD6B(19948),BCAS4(18151)-BCAS4`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`PARD6B(19948),BCAS4(18151)-BCAS4`))) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`PARD6B(19948),BCAS4(18151)-BCAS4`))) +
  theme_legend + 
  ggtitle("PARD6B(19948),BCAS4(18151)-BCAS4")  +
  scale_color_manual(values = c(
    '0' = "grey",
    '1' = "red",
    '2' = "blue",
    "3" = "green",
    "4" = "violet", 
    "na"  = "grey60"
  ))

ggsave("Fusion_PARD6B_leg.pdf", width = 10, height = 7)
ggsave("Fusion_PARD6B_noleg.pdf", width = 7, height = 7)


# Plotting NF2 fusions
FusJ <- read.csv("~/at_least_one_fusion_gene_are_coding_1249x2896.csv", check.names = FALSE, row.names = 1)
FusJ2 <- rownames_to_column(FusJ, "coordinate_ID")
finaldf <- readRDS("~/1298samples_umap2D_041923.rds")
finaldf_g <- left_join(finaldf, FusJ2, by = "coordinate_ID")
#map all NF2 fusions (high conf, recurring at least twice)
NF2fusions <- c("na")
finaldf_g$NF2fusions <- NF2fusions
finaldf_g$NF2fusions[finaldf_g$`NF2-SPATA13` >= 1] <- "NF2--SPATA13"
finaldf_g$NF2fusions[finaldf_g$`NF2-AC092595.2(36680),RN7SL252P(43770)` >= 1] <- "NF2--AC092595.2(36680),RN7SL252P(43770)"
finaldf_g$NF2fusions[finaldf_g$`NF2-CTD-2651B20.4(9455),RNU6-953P(12690)` >= 1] <- "NF2--CTD-2651B20.4(9455),RNU6-953P(12690)"
finaldf_g$NF2fusions[finaldf_g$`NF2-TCP11` >= 1] <- "NF2--TCP11"
finaldf_g$NF2fusions[finaldf_g$`NF2-TTC28` >= 1] <- "NF2--TTC28"
finaldf_g$NF2fusions[finaldf_g$`NF2-CAMK2B` >= 1] <- "NF2--CAMK2B"
finaldf_g$NF2fusions[finaldf_g$`NF2-TEX19(7298),UTS2R(3156)` >= 1] <- "NF2--TEX19(7298),UTS2R(3156)"
finaldf_g$NF2fusions[finaldf_g$`NF2-RP11-76I14.1` >= 1] <- "NF2--RP11-76I14.1"
finaldf_g$NF2fusions[finaldf_g$`NF2-NF2` >= 1] <- "NF2--NF2"
finaldf_g$NF2fusions[finaldf_g$`NF2-LARGE` >= 1] <- "NF2--LARGE"
finaldf_g$NF2fusions[finaldf_g$`NF2-CTA-85E5.10` >= 1] <- "NF2--CTA-85E5.10"
finaldf_g$NF2fusions[finaldf_g$`NF2-ASCC2` >= 1] <- "NF2--ASCC2"
finaldf_g$NF2fusions[finaldf_g$`NF2-RF00002(1304388),FRG1GP(48119)` >= 1] <- "NF2--RF00002(1304388),FRG1GP(48119)"
finaldf_g$NF2fusions[finaldf_g$`NF2-PIEZO2` >= 1] <- "NF2--PIEZO2"
finaldf_g$NF2fusions[finaldf_g$`GUCD1-NF2` >= 1] <- "GUCD1--NF2"
finaldf_g$NF2fusions[finaldf_g$`NF2-GUCD12` >= 1] <- "NF2--GUCD1"
finaldf_g$NF2fusions[finaldf_g$`DENR-NF2` >= 1] <- "DENR--NF2"

finaldf_g <- finaldf_g[order(sub("na", "", finaldf_g$NF2fusions)), ]

ggplot( finaldf_g,
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(NF2fusions),
        )
) + 
  geom_point(size=2) +
  theme_nolegend +
  ggtitle("NF2 Fusions")  +
  scale_color_manual(values = c(
    "NF2--SPATA13" = "black", 
    "NF2--AC092595.2(36680),RN7SL252P(43770)" = "red", 
    "NF2--CTD-2651B20.4(9455),RNU6-953P(12690)" = "blue", 
    "NF2--TCP11" = "sienna1", 
    "NF2--TTC28" = "green", 
    "NF2--CAMK2B" = "cyan",
    "NF2--TEX19(7298),UTS2R(3156)" = "red4",
    "NF2--RP11-76I14.1" = "orange",
    "NF2--NF2" = "steelblue1",
    "NF2--LARGE" = "hotpink4",
    "NF2--CTA-85E5.10" = "gold",
    "NF2--ASCC2" = "yellow4",
    "NF2--RF00002(1304388),FRG1GP(48119)" = "forestgreen",
    "NF2--PIEZO2" = "lightcoral",
    "GUCD1--NF2" = "lightseagreen",
    "NF2--GUCD1" = "purple",
    "DENR--NF2" = "hotpink",
    "na"  = "grey"
  ))
ggsave("NF2_fusions_RNAseq__leg.pdf", width = 7, height = 7)
ggsave("NF2_fusions_RNAseq__noleg.pdf", width = 7, height = 7)
