#title: Coloring in UMAP by main 9 clusters and subclusters. Generating Kaplan Miere plots for main clusters and subclusters.
#authors: Nayanga Thirimanne


install.packages("survival")
library(survival)
install.packages("ggsurvfit")
library(ggsurvfit)
install.packages("survminer")
library(survminer)
library(ggplot2)
library(xlsx)
library(dplyr)

sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1,]
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
#remove tumors that have progressed instead of recurred
#remove tumors that don't have data on time to recurrence
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))

time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))
timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)

# KM plots for main clusters (clusters identified using dbscan and unclustered tumors were removed)
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

timerec$cluster[timerec$coordinate_ID %in% cluster_1$coordinate_ID] <- "A"
timerec$cluster[timerec$coordinate_ID %in% cluster_2$coordinate_ID] <- "B"
timerec$cluster[timerec$coordinate_ID %in% cluster_3$coordinate_ID] <- "C"
timerec$cluster[timerec$coordinate_ID %in% cluster_4$coordinate_ID] <- "D"
timerec$cluster[timerec$coordinate_ID %in% cluster_5$coordinate_ID] <- "E"
timerec$cluster[timerec$coordinate_ID %in% cluster_6$coordinate_ID] <- "F"
timerec$cluster[timerec$coordinate_ID %in% cluster_7$coordinate_ID] <- "G"
#timerec$cluster[timerec$coordinate_ID %in% cluster_8$coordinate_ID] <- "H"
#timerec$cluster[timerec$coordinate_ID %in% cluster_9$coordinate_ID] <- "I"

timerec <- na.omit(timerec)

timerec$AB[timerec$cluster == "A"] <- "A"
timerec$AB[timerec$cluster == "B"] <- "B"
timerec$BC[timerec$cluster == "C"] <- "C"
timerec$BC[timerec$cluster == "B"] <- "B"
timerec$DE[timerec$cluster == "D"] <- "D"
timerec$DE[timerec$cluster == "E"] <- "E"
timerec$EF[timerec$cluster == "F"] <- "F"
timerec$EF[timerec$cluster == "E"] <- "E"
timerec$CF[timerec$cluster == "C"] <- "C"
timerec$CF[timerec$cluster == "F"] <- "F"
timerec$CD[timerec$cluster == "C"] <- "C"
timerec$CD[timerec$cluster == "D"] <- "D"
timerec$CE[timerec$cluster == "C"] <- "C"
timerec$CE[timerec$cluster == "E"] <- "E"
timerec$BD[timerec$cluster == "D"] <- "D"
timerec$BD[timerec$cluster == "B"] <- "B"
timerec$BE[timerec$cluster == "E"] <- "E"
timerec$BE[timerec$cluster == "B"] <- "B"
timerec$BF[timerec$cluster == "F"] <- "F"
timerec$BF[timerec$cluster == "B"] <- "B"
timerec$AG[timerec$cluster == "G"] <- "G"
timerec$AG[timerec$cluster == "A"] <- "A"
timerec$AC[timerec$cluster == "A"] <- "A"
timerec$AC[timerec$cluster == "C"] <- "C"
timerec$G[timerec$cluster == "G"] <- "G"

##
fit1 <- surv_fit(Surv(time, status) ~ AB, data = timerec)
fit2 <- surv_fit(Surv(time, status) ~ BD, data = timerec)
fit3 <- surv_fit(Surv(time, status) ~ BC, data = timerec)
fit4 <- surv_fit(Surv(time, status) ~ DE, data = timerec)
fit5 <- surv_fit(Surv(time, status) ~ EF, data = timerec)
fit6 <- surv_fit(Surv(time, status) ~ BE, data = timerec)
fit7 <- surv_fit(Surv(time, status) ~ CF, data = timerec)
fit8 <- surv_fit(Surv(time, status) ~ BF, data = timerec)
fit9 <- surv_fit(Surv(time, status) ~ G, data = timerec)
fit10 <- surv_fit(Surv(time, status) ~ AG, data = timerec)
fit11 <- surv_fit(Surv(time, status) ~ AC, data = timerec)
fit12 <- surv_fit(Surv(time, status) ~ CD, data = timerec)
fit13 <- surv_fit(Surv(time, status) ~ CE, data = timerec)

fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)

fit.list <- list(AB = fit1, BD = fit2, BC = fit3, DE = fit4, EF = fit5, BE = fit6, CF = fit7, BF = fit8, 
                 G = fit9, AG = fit10, AC = fit11, CD = fit12, CE = fit13)
surv_pvalue(fit.list)


plot_title_size = 20
legend_pt_size = 4
axis_text_size = 15
axis_title_size = 18
legend_text_size = 18
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_KM <- theme_bw()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_blank() ,
        axis.text = element_text(size = axis_text_size, face = "bold"),
        axis.title = element_text(size = axis_title_size, face = "bold")
  )
mycol <- c("Green", "Blue", "Red", "Black", "Brown", "Gold", "Cyan", "Purple")

ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = FALSE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)",ggtheme = theme_KM , palette = mycol) 
ggsave("KM_allmajorclusters.pdf", width = 7, height = 7)

####
#KM plots for subclusters in each major cluster
####

##
#subclusters in cluster A
##

sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v45_20230913.xlsx")
#sample_info <- sample_info[-1, ]
colnames(sample_info)[2] <- "sampleID_new"
sampleDeident <- read.csv("~/Deidentification/SampleID_new.csv")
sample_info_ID <- left_join(sample_info, sampleDeident, by="sampleID_new")
sample_1298 <- subset(sample_info_ID, sample_info_ID$coordinate_ID %in% colnames(log2tpm))
#remove tumors that have progressed instead of recurred
#remove tumors that don't have data on time to recurrence
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))
time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))
timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)
timerec <- na.omit(timerec)

#load samples in each subcluster 
#samples in subclusters were identified and downloaded via Oncoscape
A3Onc <- read.csv("~/A3_Onc_121923.csv")
A2Onc <- read.csv("~//A2_Onc_121923.csv")
A4Onc <- read.csv("~/A4_Onc_121923.csv")

A3Onc$Patient <- toupper(A3Onc$Patient)
A2Onc$Patient <- toupper(A2Onc$Patient)
A4Onc$Patient <- toupper(A4Onc$Patient)

A3Onc <- subset(sample_1298, sample_1298$sampleID_new %in% A3Onc$Patient)
A2Onc <- subset(sample_1298, sample_1298$sampleID_new %in% A2Onc$Patient)
A4Onc <- subset(sample_1298, sample_1298$sampleID_new %in% A4Onc$Patient)

#combine A2, A3 and A4
Asub <- c(A3Onc$coordinate_ID, A2Onc$coordinate_ID, A4Onc$coordinate_ID)
#subtract A2+A3+A4 from cluster A to identify subcluster A1
A1 <- subset(cluster_1, !(cluster_1$coordinate_ID %in% Asub))

timerec$cluster[timerec$coordinate_ID %in% A1$coordinate_ID] <- "A1"
timerec$cluster[timerec$coordinate_ID %in% A3Onc$coordinate_ID] <- "A3"
timerec$cluster[timerec$coordinate_ID %in% A2Onc$coordinate_ID] <- "A2"
timerec$cluster[timerec$coordinate_ID %in% A4Onc$coordinate_ID] <- "A4"

timerec$A12[timerec$cluster == "A1"] <- "A1"
timerec$A12[timerec$cluster == "A2"] <- "A2"
timerec$A13[timerec$cluster == "A1"] <- "A1"
timerec$A13[timerec$cluster == "A3"] <- "A3"
timerec$A23[timerec$cluster == "A2"] <- "A2"
timerec$A23[timerec$cluster == "A3"] <- "A3"
timerec$A34[timerec$cluster == "A3"] <- "A3"
timerec$A34[timerec$cluster == "A4"] <- "A4"

fit1 <- surv_fit(Surv(time, status) ~ A12, data = timerec)
fit2 <- surv_fit(Surv(time, status) ~ A13, data = timerec)
fit3 <- surv_fit(Surv(time, status) ~ A23, data = timerec)
fit4 <- surv_fit(Surv(time, status) ~ A34, data = timerec)

fit.list <- list(A13 = fit2, A12 = fit1, A23 = fit3, A34 = fit4)
surv_pvalue(fit.list)

summary(survfit2(Surv(time, status) ~ cluster, data = timerec), times = 60)
surv_pvalue(fitAll)

fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)
mycol <- c("Gold", "Blue", "Red", "Purple", "Brown", "Yellow", "Cyan", "Purple")
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = FALSE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)", ggtheme = theme_KM, palette = mycol)
ggsave("KM_subcluster_A_v3.pdf", width = 7, height = 7)

##
#subclusters in cluster C
##

sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1, ]
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
#remove tumors that have progressed instead of recurred
#remove tumors that don't have data on time to recurrence
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))
time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))
timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)
timerec <- na.omit(timerec)

C1Onc <- read.csv("~/C1Onc_36.csv")
C2Onc <- read.csv("~/C2Onc_7.csv")
C3Onc <- read.csv("~/C3Onc_23.csv")
C4Onc <- read.csv("~/C4Onc_129.csv")
C1Onc$Patient <- toupper(C1Onc$Patient)
C2Onc$Patient <- toupper(C2Onc$Patient)
C3Onc$Patient <- toupper(C3Onc$Patient)
C4Onc$Patient <- toupper(C4Onc$Patient)

timerec$cluster[timerec$coordinate_ID %in% C1Onc$Patient] <- "C1"
timerec$cluster[timerec$coordinate_ID %in% C2Onc$Patient] <- "C2"
timerec$cluster[timerec$coordinate_ID %in% C3Onc$Patient] <- "C3"
timerec$cluster[timerec$coordinate_ID %in% C4Onc$Patient] <- "C4"

fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)
mycol <- c("Black", "Red", "Blue", "green", "Brown", "Yellow", "Cyan", "Purple")

timerec$C12[timerec$cluster == "C1"] <- "C1"
timerec$C12[timerec$cluster == "C2"] <- "C2"
timerec$C34[timerec$cluster == "C3"] <- "C3"
timerec$C34[timerec$cluster == "C4"] <- "C4"
timerec$C14[timerec$cluster == "C1"] <- "C1"
timerec$C14[timerec$cluster == "C4"] <- "C4"
timerec$C13[timerec$cluster == "C1"] <- "C1"
timerec$C13[timerec$cluster == "C3"] <- "C3"

fit1 <- surv_fit(Surv(time, status) ~ C12, data = timerec)
fit2 <- surv_fit(Surv(time, status) ~ C34, data = timerec)
fit3 <- surv_fit(Surv(time, status) ~ C14, data = timerec)
fit4 <- surv_fit(Surv(time, status) ~ C13, data = timerec)
fit.list <- list(C34 = fit2, C12 = fit1, C14 = fit3, C13 = fit4)
surv_pvalue(fit.list)

ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = FALSE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)", ggtheme = theme_KM, palette = mycol)
ggsave("KM_Subcluster_C.pdf", width = 7, height = 7)

##
#subclusters in cluster F
##

sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1, ]
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
#remove tumors that have progressed instead of recurred
#remove tumors that don't have data on time to recurrence
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))
time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))
timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)
timerec <- na.omit(timerec)

F1Onc <- read.csv("~/F1Onc_19.csv")
F2Onc <- read.csv("~/F2Onc_39.csv")
F1Onc$Patient <- toupper(F1Onc$Patient)
F2Onc$Patient <- toupper(F2Onc$Patient)

timerec$cluster[timerec$coordinate_ID %in% F1Onc$Patient] <- "F1"
timerec$cluster[timerec$coordinate_ID %in% F2Onc$Patient] <- "F2"

fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)
mycol <- c("Red", "Blue")
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = FALSE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)", ggtheme = theme_KM, palette = mycol)
ggsave("KM_Subcluster_F.pdf", width = 7, height = 7)

##
#subclusters in cluster B
##

sample_info <- readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1, ]
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
#remove tumors that have progressed instead of recurred
#remove tumors that don't have data on time to recurrence
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))
time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))
timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)
timerec <- na.omit(timerec)

B1Onc <- read.csv("~/B1Onc_196.csv")
B2Onc <- read.csv("~/B2Onc_75.csv")
B1Onc$Patient <- toupper(B1Onc$Patient)
B2Onc$Patient <- toupper(B2Onc$Patient)
timerec$cluster[timerec$coordinate_ID %in% B1Onc$Patient] <- "B1"
timerec$cluster[timerec$coordinate_ID %in% B2Onc$Patient] <- "B2"
fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)
mycol <- c("Blue", "Red")
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)", ggtheme = theme_KM, palette = mycol)
ggsave("KM_Subcluster_B.pdf", width = 7, height = 7)


####
#color subclusters in UMAP 
####

## subclusters in cluster A
finaldf_A <- finaldf
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% A3Onc$coordinate_ID] <- "A3"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% A2Onc$coordinate_ID] <- "A2"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% A1$coordinate_ID] <- "A1"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% A4Onc$coordinate_ID] <- "A4"
finaldf_A$Afinal[is.na(finaldf_A$Afinal)] <- "na"

ggplot(
  finaldf_A,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = Afinal
  )
) + 
  geom_point(size=1.5) +
  theme_legend +
  scale_color_manual(values = c("A1" = "gold", 
                                "A2" = "blue", 
                                "A3" = "red",
                                "A4" = "purple",
                                "na" = "grey"))
ggsave("Subcluster_A_UMAP_2.pdf", width = 7, height = 7)
ggsave("Subcluster_A_UMAP_leg_2.pdf", width = 10, height = 7)

## subclusters in cluster C
finaldf_A <- finaldf
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% C1Onc$Patient] <- "C1"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% C2Onc$Patient] <- "C2"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% C3Onc$Patient] <- "C3"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% C4Onc$Patient] <- "C4"
finaldf_A$Afinal[is.na(finaldf_A$Afinal)] <- "na"

ggplot(
  finaldf_A,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = Afinal
  )
) + 
  geom_point(size=1.5) +
  theme_legend +
  scale_color_manual(values = c("C1" = "black",
                                "C2" = "red", 
                                "C3" = "green",
                                "C4" = "blue", 
                                "na" = "grey"))
ggsave("Subcluster_C_UMAP.pdf", width = 7, height = 7)
ggsave("Subcluster_C_UMAP_leg.pdf", width = 10, height = 7)

## subclusters in cluster B
finaldf_A <- finaldf
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% B1Onc$Patient] <- "B1"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% B2Onc$Patient] <- "B2"
finaldf_A$Afinal[is.na(finaldf_A$Afinal)] <- "na"
ggplot(
  finaldf_A,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = Afinal
  )
) + 
  geom_point(size=1.5) +
  theme_legend +
  scale_color_manual(values = c("B1" = "blue", 
                                "B2" = "red", 
                                "na" = "grey"))
ggsave("Subcluster_B_UMAP.pdf", width = 7, height = 7)
ggsave("Subcluster_B_UMAP_leg.pdf", width = 10, height = 7)


## subclusters in cluster F 
finaldf_A <- finaldf
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% F1Onc$Patient] <- "F1"
finaldf_A$Afinal[finaldf_A$coordinate_ID %in% F2Onc$Patient] <- "F2"
finaldf_A$Afinal[is.na(finaldf_A$Afinal)] <- "na"
ggplot(
  finaldf_A,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = Afinal
  )
) + 
  geom_point(size=1.5) +
  theme_legend +
  scale_color_manual(values = c("F1" = "red", 
                                "F2" = "blue", 
                                "na" = "grey"))
ggsave("Subcluster_F_UMAP.pdf", width = 7, height = 7)
ggsave("Subcluster_F_UMAP_leg.pdf", width = 10, height = 7)
