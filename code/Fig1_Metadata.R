#title: Generate UMAP and color in UMAP by metadata
#authors: Nayanga Thirimanne

library(ggplot2)
install.packages("umap")
library(umap)
library(tibble)
library(readxl)
library(dplyr)

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
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_blank())

############################
# Calculate UMAP coordinates
############################

log2_tpm_adjusted <- #load batch corrected, log2tpm counts
umap_out <- umap(t(log2_tpm_adjusted), random_state = 123, min_dist = 0.1, metric = 'cosine')
umap_2d <- data.frame(umap_out$layout) %>%
  tibble::rownames_to_column("samples")
umap_2d$samples <- gsub("^X", "", umap_2d$samples)
colnames(umap_2d) = c("coordinate_ID","UMAP1_2D", "UMAP2_2D")
# Merge with metadata
sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip = 5)
sample_info <- sample_info[-1, ]
sample_info_2 <- sample_info[match(umap_2d$coordinate_ID, sample_info$coordinate_ID), ] 
finaldf = left_join(umap_2d, sample_info_2,by="coordinate_ID" )
finaldf <- write.csv(finaldf, "~/1298samples_umap2D_041923.csv", check.names = FALSE)

##########################
#Reproduce figure 1 panels
##########################

# WHO Grade
ggplot( finaldf %>% 
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
    'WHO 1' = "orange",
    'WHO 2' = "forestgreen",
    'WHO 3' = "red",
    "na" = "grey",
    "human YAP1fus" = "blue"
  ))

# Recurrent tumors
ggplot( finaldf %>% 
          arrange(factor(RECURRENT_TUMOR)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(RECURRENT_TUMOR),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("Recurrent Tumors")  +
  scale_color_manual(values = c(
    'yes' = "blue",
    'no' = "orange",
    "na" = "grey"
  ))

# CHR 22 Loss
ggplot( finaldf %>% 
          arrange(factor(CHR22LOSS)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(CHR22LOSS),
        )
) + 
  geom_point(size=1.5) +
  theme_nolegend + 
  ggtitle("Chr 22 Loss") +
  scale_color_manual(values = c('yes' = "blue",
                                'no' = "orange",
                                "na" = "grey" 
  ))

# NF2 mutations
ggplot( finaldf %>% 
          arrange(factor(NF2)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(NF2),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("NF2 Mutations") +
  scale_color_manual(values = c('yes' = "blue",
                                'no' = "orange",
                                "na" = "grey"))


# NF2 Fusions
ggplot( finaldf %>% 
          arrange(factor(NF2_FUSION)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(NF2_FUSION),
        )
) + 
  geom_point(size=1.5) +
  theme_legend +
  ggtitle("NF2 Fusions") +
  scale_color_manual(values = c('yes' = "blue",
                                'no' = "orange",
                                "na" = "grey"))

# NF2 Mutations + NF2 Fusions + Chr22 Loss
ggplot( finaldf %>% 
          arrange(factor(NF2_ISSUE)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(NF2_ISSUE),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("NF2 Mut + NF2 Fusions + Chr 22 Loss ")  +
  scale_color_manual(values = c('yes' = "blue",
                                'no' = "orange",
                                "na" = "grey"))

# Gender
MALE <- subset(finaldf, finaldf$GENDER == "male")
ggplot( finaldf %>% 
          arrange(desc(factor(GENDER))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(GENDER),
        )
) + 
  geom_point(size=1.5) +
  geom_point(data=MALE,
             aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(GENDER)))+
  theme_legend + 
  ggtitle("Gender") +
  scale_color_manual(values = c("male" = "blue",
                                "female" = "violet",
                                "na" = "grey"))

# TRAF7, AKT1, KLF4 and SMO mutations
ggplot( finaldf %>% 
          arrange(factor(TRAF7_AKT1_KLF4_SMO)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(TRAF7_AKT1_KLF4_SMO),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("TRAF7 AKT1 KLF4 SMO") +
  scale_color_manual(values = c('yes' = "blue",
                                'no' = "orange",
                                "na" = "grey"))

# Tumors with more than one of TRAF7, KLF4, AKT1, SMO mutations
finaldf$nonNF2[finaldf$TRAF7 == "na"] <- "na"
finaldf$nonNF2[finaldf$KLF4 == "na"] <- "na"
finaldf$nonNF2[finaldf$SMO == "na"] <- "na"
finaldf$nonNF2[finaldf$AKT1 == "na"] <- "na"
finaldf$nonNF2[finaldf$TRAF7 == "no"] <- "na"
finaldf$nonNF2[finaldf$KLF4 == "no"] <- "na"
finaldf$nonNF2[finaldf$SMO == "no"] <- "na"
finaldf$nonNF2[finaldf$AKT1 == "no"] <- "na"
finaldf$nonNF2[finaldf$TRAF7 == "yes"] <- "TRAF7"
finaldf$nonNF2[finaldf$KLF4 == "yes"] <- "KLF4"
finaldf$nonNF2[finaldf$SMO == "yes"] <- "SMO"
finaldf$nonNF2[finaldf$AKT1 == "yes"] <- "AKT1"
finaldf$nonNF2[finaldf$TRAF7 == "yes" & finaldf$KLF4=="yes"] <- "TRAF7, KLF4"
finaldf$nonNF2[finaldf$TRAF7 == "yes" & finaldf$AKT1=="yes"] <- "TRAF7, AKT1"
finaldf$nonNF2[finaldf$TRAF7 == "yes" & finaldf$SMO=="yes"] <- "TRAF7, SMO"
finaldf$nonNF2[finaldf$TRAF7 == "yes" & finaldf$KLF4=="yes" & finaldf$SMO=="yes"] <- "TRAF7, KLF4, SMO"
finaldf$nonNF2[finaldf$TRAF7 == "yes" & finaldf$NF2=="yes"] <- "TRAF7, NF2"

ggplot( finaldf %>% 
          arrange(factor(nonNF2)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(nonNF2),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("TRAF7 AKT1 KLF4 SMO") +
  scale_color_manual(values = c("TRAF7" = "blue",
                                "KLF4" = "red",
                                "AKT1" = "green",
                                "SMO" = "brown",
                                "TRAF7, KLF4" = "purple",
                                "TRAF7, AKT1" = "deepskyblue",
                                "TRAF7, SMO" = "black",
                                "TRAF7, KLF4, SMO" = "violet",
                                "TRAF7, NF2" = "cyan",
                                "no" = "gold",
                                "na" = "grey"))


# SMO and TRAF7 mutations
finaldf$SMO_TRAF7 <- "na"
finaldf$SMO_TRAF7[finaldf$SMO == "na"] <- "na"
finaldf$SMO_TRAF7[finaldf$SMO == "no"] <- "na"
finaldf$SMO_TRAF7[finaldf$TRAF7 == "na"] <- "na"
finaldf$SMO_TRAF7[finaldf$TRAF7 == "no"] <- "na"
finaldf$SMO_TRAF7[finaldf$SMO == "yes"] <- "SMO"
finaldf$SMO_TRAF7[finaldf$TRAF7 == "yes" & finaldf$SMO == "yes"] <- "TRAF7,SMO"

SMO <- subset(finaldf, finaldf$SMO_TRAF7 == "SMO")
SMOTRAF7 <- subset(finaldf, finaldf$SMO_TRAF7 == "TRAF7,SMO")

ggplot( finaldf %>% 
          arrange(factor(SMO_TRAF7)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(SMO_TRAF7),
        )
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("SMO mutations") +
  scale_color_manual(values = c(
    "SMO" = "brown",
    "TRAF7,SMO" = "black",
    "na" = "grey"))


#KLF4 and AKT1 mutations
finaldf$KLF4AKT1[finaldf$TRAF7 == "na"] <- "na"
finaldf$KLF4AKT1[finaldf$KLF4 == "na"] <- "na"
finaldf$KLF4AKT1[finaldf$SMO == "na"] <- "na"
finaldf$KLF4AKT1[finaldf$AKT1 == "na"] <- "na"
finaldf$KLF4AKT1[finaldf$TRAF7 == "no"] <- "na"
finaldf$KLF4AKT1[finaldf$KLF4 == "no"] <- "na"
finaldf$KLF4AKT1[finaldf$SMO == "no"] <- "na"
finaldf$KLF4AKT1[finaldf$AKT1 == "no"] <- "na"
finaldf$KLF4AKT1[finaldf$KLF4 == "yes"] <- "KLF4"
finaldf$KLF4AKT1[finaldf$AKT1 == "yes"] <- "AKT1"
KLF4 <- subset(finaldf, finaldf$KLF4AKT1 == "KLF4")
AKT1 <- subset(finaldf, finaldf$KLF4AKT1 == "AKT1")

ggplot( finaldf %>% 
          arrange(factor(KLF4AKT1)),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(KLF4AKT1),
        )
) + 
  geom_point(size=1.5) +
  geom_point(data=KLF4, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(KLF4AKT1))) +
  geom_point(data=AKT1, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(KLF4AKT1))) +
  theme_nolegend + 
  ggtitle("KLF4 and AKT1 mutations") +
  scale_color_manual(values = c(
    "KLF4" = "red",
    "AKT1" = "blue",
    "na" = "grey"))


# YAP1 fusions
MAML <- subset(finaldf, finaldf$YAP1_Fusion == "YAP1-MAML2")
LMO <- subset(finaldf, finaldf$YAP1_Fusion == "YAP1-LMO1")
ggplot( finaldf %>% 
          arrange(desc(factor(YAP1_Fusion))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(YAP1_Fusion),
        )
) + 
  geom_point(size=1.5) +
  geom_point(data=MAML,
             aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(YAP1_Fusion)))+
  geom_point(data=LMO,
             aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(YAP1_Fusion)))+
  theme_legend + 
  ggtitle("YAP1 Fusions") +
  scale_color_manual(values = c("YAP1-MAML2" = "red",
                                "YAP1-LMO1" = "blue",
                                "na" = "grey"))


# Age categories 
finaldf$AGE <- as.numeric(finaldf$AGE)
finaldf$AGE_CATEGORY <- cut(finaldf$AGE, breaks = c(-1, 0, 9, 18, 30, 60, 90, 98),
                            labels = c("1--9", "1--9", "10--18", "19--30","31--60", "61--90", "91--98"))

finaldf$AGE_CATEGORY <- as.character(finaldf$AGE_CATEGORY)
finaldf$AGE_CATEGORY[is.na(finaldf$AGE_CATEGORY)] <- "na"
finaldf$AGE_CATEGORY <- as.factor(finaldf$AGE_CATEGORY)

One <- subset(finaldf, finaldf$AGE_CATEGORY == "1--9")
Two <- subset(finaldf, finaldf$AGE_CATEGORY == "10--18")
Three <- subset(finaldf, finaldf$AGE_CATEGORY == "19--30")
Four <- subset(finaldf, finaldf$AGE_CATEGORY == "31--60")
Five <- subset(finaldf, finaldf$AGE_CATEGORY == "61--90")
Six <- subset(finaldf, finaldf$AGE_CATEGORY == "91--98")
Seven <- subset(finaldf, finaldf$AGE_CATEGORY == "na")

ggplot(finaldf %>% 
         arrange(factor(AGE_CATEGORY)), aes(x = UMAP1_2D, y = UMAP2_2D)) + 
  geom_point(size=1.5) +
  theme_nolegend +
  ggtitle("Age") + 
  geom_point(data = Seven, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = Four, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = Five, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = Six, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = Three, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = Two, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  geom_point(data = One, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(AGE_CATEGORY))) +
  scale_color_manual(values = c(
    '1--9' = " blue",
    '10--18' = "firebrick1 ",
    "19--30" = "sandybrown",
    "31--60" = " khaki3",
    "61--90" = " lavenderblush3", 
    "91--98" = "paleturquoise4 ", 
    "na"  = "grey39"
  ))


# Bar graph showing Age Categories in 2 regions (Region 2 and rest of the UMAP) (Supplemental Figure)
Reg1 <- read.csv("~/Analysis30_Age_Gender/Region_younger.csv") #downloaded samples in Region 2 from Oncoscape
Reg1$Patient <- toupper(Reg1$Patient)

# Assign Region labels based on coordinate_ID
finaldf$Region[finaldf$coordinate_ID %in% Reg1$Patient] <- "Region 2"
finaldf$Region[!(finaldf$coordinate_ID %in% Reg1$Patient)] <- "All-Region2"

finaldf <- finaldf[!is.na(finaldf$Region), ]

finaldf$AGE_CATEGORY <- cut(finaldf$AGE, breaks = c(-1, 0, 9, 18, 30, 60, 90, 98),
                            labels = c("1--9", "1--9", "10--18", "19--30","31--60", "61--90", "91--98"))
finaldf$AGE_CATEGORY <- as.character(finaldf$AGE_CATEGORY)
finaldf$AGE_CATEGORY[is.na(finaldf$AGE_CATEGORY)] <- "na"
finaldf$AGE_CATEGORY <- as.factor(finaldf$AGE_CATEGORY)

#Group and summarize data
BG_tbl <- finaldf %>% group_by(Region, AGE_CATEGORY) %>% 
  summarise(total_count=n(),
            .groups = "drop") %>% 
  as.data.frame()

#Calculate counts in Region X and All-X
Region_X_count <- sum(BG_tbl[which(BG_tbl$Region == "Region 2"), 3])
Region_AllX_count <- sum(BG_tbl[which(BG_tbl$Region == "All-Region2"), 3])

#calculate percentages for each category
BG_tbl_2 <- BG_tbl[-c(7,14), ] #remove na
BG_tbl_3 <- BG_tbl_2 %>% filter(Region != "All-Region2") %>% mutate(Percentage = total_count/sum(total_count)*100)
BG_tbl_4 <- BG_tbl_2 %>% filter(Region != "Region 2") %>% mutate(Percentage = total_count/sum(total_count)*100)
BG_tbl_5 <- rbind(BG_tbl_3, BG_tbl_4)

#Generate bar plot
ggplot(BG_tbl_5, aes(x=factor(Region, level=c("Region 2", "All-Region2")), y=Percentage))+
  geom_bar(aes(fill=AGE_CATEGORY), stat="identity", position= position_dodge(preserve = "single"), color="black") +
  theme_grid_nolegend +
  ggtitle("AGE")+
  scale_fill_manual(values = c('1--9' = " blue",
                               '10--18' = "firebrick1 ",
                               "19--30" = "sandybrown",
                               "31--60" = " khaki3",
                               "61--90" = " lavenderblush3", 
                               "91--98" = "paleturquoise4 ", 
                               "na"  = "grey39"
  )) +
  xlab("Regions")



#Gender in 2 regions (Region 1 and rest of the UMAP) 
finaldf <- read.csv("~/1298samples_umap2D_041923.csv")
Reg1 <- read.csv("~/Analysis30_Age_Gender/Region_male.csv")
Reg1$Patient <- toupper(Reg1$Patient)

finaldf$Region[finaldf$coordinate_ID %in% Reg1$Patient] <- "Region 1"
finaldf$Region[!(finaldf$coordinate_ID %in% Reg1$Patient)] <- "All-Region1"
finaldf <- finaldf[!is.na(finaldf$Region), ]

#Group and summarize data
BG_tbl <- finaldf %>% group_by(Region, GENDER) %>% 
  summarise(total_count=n(),
            .groups = "drop") %>% 
  as.data.frame()

#Calculate counts for Region X and All-X
Region_X_count <- sum(BG_tbl[which(BG_tbl$Region == "Region X"), 3])
Region_AllX_count <- sum(BG_tbl[which(BG_tbl$Region == "All-X"), 3])

# Create table without rows 3 and 6 (remove na counts)
BG_tbl_2 <- BG_tbl[-c(3,6), ]
perc <- c(69,31,39,61) #manually calculated percentages
BG_tbl_2$percentage <- perc

#Create a bar plot (Supplemental Figure)
ggplot(BG_tbl_2, aes(x=factor(Region, level=c("Region 1", "All-Region1")), y=percentage))+
  geom_bar(aes(fill=GENDER), stat="identity", position= position_dodge(preserve = "single"), color="black") +
  theme_grid_legend +
  ggtitle("Gender")+
  scale_fill_manual(values = c('male' = "blue",
                               'female' = "violet",
                               "na"  = "grey"
  )) +
  xlab("Regions")

## FFPE and Frozen tissue
ggplot( finaldf %>% 
  arrange(desc(factor(FFPE_or_Frozen))),
aes(
  x = UMAP1_2D, 
  y = UMAP2_2D,
  color = factor(FFPE_or_Frozen),
)
) + 
  geom_point(size=1.5) +
  theme_legend + 
  ggtitle("FFPE and frozen samples") +
  scale_color_manual(values = c("FFPE" = "blue",
                                "frozen" = "grey80",
                                "na" = "grey"))
ggsave("FFPE_frozen_leg.pdf", width = 10, height = 7)
ggsave("FFPE_frozen_noleg.pdf", width = 7, height = 7)


## Months to recurrence
sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1,]
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))
#sample_timetorec <- sample_1298
sample_timetorec$time_for_KM_final <- as.numeric(as.character(sample_timetorec$time_for_KM_final))

sample_timetorec$time_for_KM_final_category <- cut(sample_timetorec$time_for_KM_final, breaks = c(-1, 0, 24, 60, 120, 315),
                                                   labels = c("0-24", "0-24", "25-60", "61-120", "120<"))

sample_info_sub <- sample_timetorec[, c("coordinate_ID", "time_for_KM_final_category", "has_the_tumor_reccur", "time_for_KM_final")]
sample_info_sub_2 <- subset(sample_info_sub, sample_info_sub$has_the_tumor_reccur == "yes")

fianldf <- left_join(finaldf, sample_info_sub_2, by="coordinate_ID")

fianldf$time_for_KM_final_category <- as.character(fianldf$time_for_KM_final_category)
fianldf$time_for_KM_final_category[is.na(fianldf$time_for_KM_final_category)] <- "na"
fianldf$time_for_KM_final_category <- as.factor(fianldf$time_for_KM_final_category)

One <- subset(fianldf, fianldf$time_for_KM_final_category == "0-24" & fianldf$has_the_tumor_reccur == "yes")
Two <- subset(fianldf, fianldf$time_for_KM_final_category == "25-60" & fianldf$has_the_tumor_reccur == "yes")
Three <- subset(fianldf, fianldf$time_for_KM_final_category == "61-120" & fianldf$has_the_tumor_reccur == "yes")
Four <- subset(fianldf, fianldf$time_for_KM_final_category == "120<" & fianldf$has_the_tumor_reccur == "yes")
Five <- subset(fianldf, fianldf$time_for_KM_final_category == "na")
ggplot( fianldf %>% 
          arrange(desc(factor(time_for_KM_final_category))),
        aes(
          x = UMAP1_2D, 
          y = UMAP2_2D,
          color = factor(time_for_KM_final_category),
        )
) + 
  geom_point( size=1.5) +
  theme_legend +
  ggtitle("MONTHS TO RECURRENCE") +
  scale_color_manual(values = c( 
    "0-24" = "red",
    "25-60" = "gold",
    "61-120" = "green4",
    "120<" = "royalblue1", 
    "na" = "grey"))
ggsave("MonthsToRec_leg.pdf", width = 8, height = 7)
ggsave("MonthsToRec_noleg.pdf", width = 7, height = 7)
