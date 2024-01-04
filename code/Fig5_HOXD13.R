#title: Comparing tumors with high and low HOXD13 levels in cluster A
#authors: Nayanga Thirimanne

library(dplyr)
install.packages("survival")
library(survival)
install.packages("ggsurvfit")
library(ggsurvfit)
install.packages("survminer")
library(survminer)

#load normalized batch corrected count file
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")
#load dbscan cluster IDs
clusterID <- readRDS("~/dbscan_1298_cl9_uncl48.rds")
#remove unclustered samples
unclust <- clusterID[clusterID$cluster == "0",]
clusterID <- subset(clusterID, !(clusterID$coordinate_ID %in% unclust$coordinate_ID))

#select the gene of interest and convert it to a dataframe
HOXp1 <- t(log2tpm) 
HOXp <- HOXp1[, c("HOXD13")] %>% as.data.frame()
colnames(HOXp)[1] <- "HOXD13"
HOXp <- rownames_to_column(HOXp, var="coordinate_ID")

HOXpp <- left_join(HOXp, clusterID, by="coordinate_ID")

#subset cluster A samples
HOXD_A <- subset(HOXpp, HOXpp$cluster == 1)

#box plot of HOXD13 expression in cluster A 
p <- ggplot(HOXD_A, aes(x="", y=HOXD13)) +
  geom_boxplot()
p + geom_jitter() + ggtitle("HOXD13 in Cluster A")


QH <- quantile(HOXD_A$HOXD13, 0.75)
QL <- quantile(HOXD_A$HOXD13, 0.25)

HOXD13_A_H <- HOXD_A[HOXD_A$HOXD13 > QH, ]
HOXD13_A_L <- HOXD_A[HOXD_A$HOXD13 < QL, ]


#KM plots for high vs low HOXD13 samples cluster A
#read metadata file to extract time to recurrence information
sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip = 5)
sample_info <- sample_info[-1, ]
sample_1298 <- subset(sample_info, sample_info$coordinate_ID %in% colnames(log2tpm))
sample_timetorec <- subset(sample_1298, !(sample_1298$has_the_tumor_progressed == "yes")) %>% subset(., !(.$time_for_KM_final == "na"))

time_rec <- sample_timetorec[, c("coordinate_ID","time_for_KM_final", "has_the_tumor_reccur")]
colnames(time_rec) <- c("coordinate_ID", "time", "status")
#remove samples with status=na
timerec <- subset(time_rec, !(time_rec$status == "na")) %>% subset(., !(.$time == "na"))

timerec$status[timerec$status == "yes"] <- "1"
timerec$status[timerec$status == "no"] <- "0"
timerec$time <- as.numeric( timerec$time)
timerec$status <- as.numeric( timerec$status)

timerec$cluster[timerec$coordinate_ID %in% HOXD13_A_H$coordinate_ID] <- "A_High"
timerec$cluster[timerec$coordinate_ID %in% HOXD13_A_L$coordinate_ID] <- "A_low"

timerec$A12[timerec$cluster == "A_High"] <- "A_High"
timerec$A34[timerec$cluster == "A_low"] <- "A_low"

fit1 <- surv_fit(Surv(time, status) ~ A12, data = timerec)
fit2 <- surv_fit(Surv(time, status) ~ A34, data = timerec)

#timerec <- na.omit(timerec)
fitAll <- surv_fit(Surv(time, status) ~ cluster, data = timerec)
summary(survfit2(Surv(time, status) ~ cluster, data = timerec), times = 60)
surv_pvalue(fitAll)

fit.list <- list(A12 = fit1, A34 = fit2)
surv_pvalue(fit.list)

mycol <- c("Red", "Gold")
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)", ggtheme = theme_KM, palette = mycol)

ggsave("ClusterA_high_low_HOXD13_KM.pdf", height = 3, width=4)
