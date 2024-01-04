#title: Copy Number Alterations using CaSpER on bulk RNA-Seq
#authors: Nayanga Thirimanne

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")

library(sva)
library(CaSpER)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(edgeR)
install.packages("survival")
library(survival)
install.packages("ggsurvfit")
library(ggsurvfit)
install.packages("survminer")
library(survminer)

BiocManager::install("BiocGenerics")
library(BiocGenerics)
library(IRanges)

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



####
#Run CaSpER 
####

# Here we are using hippocampus and frontol cortex normal tumors as controls
# For normal samples RNA-Seq data was downloaded from GTEX and processed using the same pipelline

#normalize and batch correct meningioma datasets along with GTEX samples
df1 <- readRDS("~/Meningioma_dataset1_538/GSE139652_10/GSE139652_10samples_proteincodingGencode_rawcounts.rds")
df2 <- readRDS("~/Meningioma_dataset1_538/GSE136661_160/GSE136661_159samples_proteincodingGencode_rawcounts.rds")
df3 <- readRDS("~/Meningioma_dataset1_538/Holland_cohort_1_101/Hollandcohort1_93samples_proteincodingGencode_rawcounts.rds")
df4 <- readRDS("~/Meningioma_dataset1_538/Holland_cohort_2_108/Hollandcohort2_108samples_proteincodingGencode_rawcounts.rds")
df5 <- readRDS("~/Meningioma_dataset1_538/Felix_42/Felix_44samples_proteincodingGencode_rawcounts.rds")
df6 <- readRDS("~/Meningioma_dataset1_538/Toronto_124/Toronto_123samples_proteincodingGencode_rawcounts.rds")
df7 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE101638/GSE101638_42samples_proteincodingGencode_rawcounts.rds")
df8 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE183653/GSE183653_185samples_proteincodingGencode_rawcounts.rds")
df9 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE151921/GSE151921_13samples_proteincodingGencode_rawcounts.rds")
df10 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE85133/GSE85133_19samples_proteincodingGencode_rawcounts.rds")
df11 <- readRDS("~/Meningioma_RNASeq_dataset2/PRJNA705586/PRJNA705586_70samples_proteincodingGencode_rawcounts.rds")
df12 <- readRDS("~/Meningioma_RNASeq_dataset2/CAVATICA_CHOP_analysis/CAVATICA_25samples_proteincodingGencode_rawcounts.rds")
df13 <- readRDS("~/Meningioma_RNASeq_dataset2/HollandUW_3/HollandUW_3_78samples_wodups_rawcounts_proteincoGencode_013123.rds")
df14 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE212666_HKU/GSE212666_301samples_proteincodingGencode_rawcounts.rds")
df15 <- readRDS("~/Meningioma_RNASeq_dataset2/CBTN_youngadults/CBTN_28samples_proteincodingGencode_rawcounts.rds")

# GTEX frontal cortex =30 and GTEX hippocampus =8
df16 <- readRDS("~/GTEX/GTEX_30frontalcor_proteincodingGencode_rawcounts.rds")
df17 <- readRDS("~/GTEX/GTEX_8hippocamp_proteincodingGencode_rawcounts.rds")

DF <- cbind(df1, df2, df3, df4, df5, df6, 
            df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17)
DF <- as.data.frame(DF)

batch <- c(rep(1, ncol(df1)), 
           rep(2, ncol(df2)),  
           rep(3, ncol(df3)), 
           rep(4, ncol(df4)),
           rep(5, ncol(df5)),
           rep(6, ncol(df6)),
           rep(7, ncol(df7)),
           rep(8, ncol(df8)),
           rep(9, ncol(df9)),
           rep(10, ncol(df10)),
           rep(11, ncol(df11)), 
           rep(12, ncol(df12)),
           rep(13, ncol(df13)), 
           rep(14, ncol(df14)), 
           rep(15, ncol(df15)),
           rep(16, ncol(df16)),
           rep(17, ncol(df17))
)

#batch correction
adjusted_CNV <-  ComBat_seq(as.matrix(DF), batch=batch, group=NULL)

#normalize read counts 
gtf <- rtracklayer::import("~/Reference_Genomes/hg38_release39_053122/gencode.v39.primary_assembly.annotation.gtf")
hg38gtf <- as.data.frame(gtf)
hg38gtf_sub <- subset(hg38gtf, gene_name %in% rownames(adjusted_CNV)) 
hg38gtf_sub2 <- dplyr::filter(hg38gtf_sub, hg38gtf_sub$gene_type == "protein_coding" & hg38gtf_sub$type == "gene")
hg38gtf_sub2 <- hg38gtf_sub2 %>% arrange(gene_name) %>%
  dplyr::filter(duplicated(gene_name) == FALSE)
gene_width = hg38gtf_sub2$width
gene_name = rownames(adjusted_CNV)
y <- DGEList(counts = adjusted_CNV)
rpkm = rpkm(y, gene.length = gene_width)
rownames(rpkm) = gene_name
tpm = apply(rpkm, 2, function(x){
  (x/sum(x))*10^6
})
log2_tpm_ar =log2(tpm+1)

#### 
df50 <- log2_tpm_ar
#drop samples that doesn't have CaSpER calls
drop <- c("SRR3995987", "SRR13810536", "H486", "1288", "SRR1319539", "SRR1342045", "SRR1388305", "SRR1402900", "SRR1433796")
df50.1 <- df50[, !(colnames(df50) %in% drop)]
cnvtestdata <- as.matrix(df50.1)

#Cytoband information
cytoband <- read.delim("~/cytoBand_hg38_ucsc.txt", header = FALSE)
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

#Centromere information
centromere <- read.delim("~/centromere_hg38.txt", header = FALSE)

#control sample ids
controls <- read.csv("~/CaSpER_test/control.sample.id.csv", header = TRUE)
control.sample.ids <- controls$controls_GTEX

#Annotation data (positions of each gene along each chromosome)
annotation <- CaSpER::generateAnnotation(id_type="hgnc_symbol", genes=rownames(df50.1), ishg19 = F, centromere, host="useast.ensembl.org")

#subset data to only have genes as same as in annotation file
cnvtestdata <- cnvtestdata[match( annotation$Gene,rownames(cnvtestdata)), ]

#Reading BAFExtract output files
#baf files of all meningioma samples were in the compile_menin folder
loh <- CaSpER::readBAFExtractOutput(path="~/CaSpER_Meningioma/compile_menin/", sequencing.type = "bulk", suffix = "baf")


#generate dataframe of loh.name.mapping
files <- list.files(pattern = ".baf")
loh.name <- c(files)
sample.name <- c(files)
loh.name.mapping <- data.frame(loh.name, sample.name)


#creating CaSpER object
object <- CaSpER::CreateCasperObject(raw.data = cnvtestdata, loh.name.mapping = loh.name.mapping, sequencing.type = "bulk", cnv.scale =3, loh.scale = 3, matrix.type = "normalized",
                                     expr.cutoff = 0.1, annotation = annotation, method = "iterative", loh=loh, filter = "median", log.transformed = TRUE,
                                     control.sample.ids = control.sample.ids, cytoband = cytoband, genomeVersion = "hg38")

#pairwise comparison of scales from BAF and expression signals
final.objects <- CaSpER::runCaSpER(object, removeCentromere = T, cytoband = cytoband, method = "iterative")

#Large scale CNV summarization
finalChrMat <- CaSpER::extractLargeScaleEvents(final.objects, thr = 0.75)

#Visualization 
obj <- final.objects[[9]]
CaSpER::plotLargeScaleEvent2(chrMat = finalChrMat, fileName = "large.scale.events.1328samples_ChrMat_050223.pdf")
dev.off()

#Visualize only GTEX samples
finalChrMat_GTEX <- finalChrMat[rownames(finalChrMat) %in% controls$controls_GTEX, ]
CaSpER::plotLargeScaleEvent2(chrMat = finalChrMat_GTEX, fileName = "large.scale.events.1328samples_GTEX_ChrMat_050223.pdf")

# color in UMAP with CNV calls
CNVcalls <- as.data.frame(finalChrMat)
# CNV burden
CNVcalls$del_calls <- rowSums(finalChrMat == "-1")
CNVcalls$amp_calls <- rowSums(finalChrMat == "1")

CNVcalls <- rownames_to_column(CNVcalls, "coordinate_ID")
colnames(CNVcalls)[1] <- "coordinate_ID"
finaldf <- readRDS("~/1298_UMAP2D.rds") #load file with UMAP 2D coordinates
finaldf_CNV = left_join(finaldf, CNVcalls,by="coordinate_ID" )


####
#Color UMAP by CNAs
####

#CNA burden - deletions
cohort4 <- subset(finaldf_CNV, finaldf_CNV$del_calls >= 0)
ggplot(
  finaldf_CNV,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D
  )
) + 
  geom_point(size=1.5, colour = "grey80") +
  geom_point(data = cohort4,
             aes(x=UMAP1_2D, y=UMAP2_2D, color=del_calls),  size = 1.5) +
  theme_legend + 
  scale_colour_gradientn(colours = c(low="lightyellow", high="red3")) +
  ggtitle("CNA burden: deletions")
ggsave("CNVburden_del_leg.pdf", width = 10, height = 7)
ggsave("CNVburden_del_noleg.pdf", width = 7, height = 7)


# Chr 10q
finaldf_CNV$`10q`[is.na(finaldf_CNV$`10q`)] <- "na"
FF <- subset(finaldf_CNV, finaldf_CNV$`10q` == 0) 
FF2 <- subset(finaldf_CNV, finaldf_CNV$`10q` == 1) 
FF3 <- subset(finaldf_CNV, finaldf_CNV$`10q` == -1) 
ggplot(
  finaldf_CNV %>% arrange(factor(`10q`)),
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(`10q`)
  )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`10q`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`10q`))) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`10q`))) +
  theme_legend +
  ggtitle("Chr 10q") +
  scale_color_manual(values = c("-1" = "blue",
                                "1" = "red",
                                "0" = "grey", 
                                "na" = "grey"))

ggsave("CNV_10q_leg.pdf", width = 10, height = 7)
ggsave("CNV_10q_noleg.pdf", width = 7, height = 7)


# Chr 1p
finaldf_CNV$`1p`[is.na(finaldf_CNV$`1p`)] <- "na"
FF <- subset(finaldf_CNV, finaldf_CNV$`1p` == 0) 
FF2 <- subset(finaldf_CNV, finaldf_CNV$`1p` == 1) 
FF3 <- subset(finaldf_CNV, finaldf_CNV$`1p` == -1) 
ggplot(
  finaldf_CNV %>% arrange(factor(`1p`)),
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(`1p`)
  )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`1p`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`1p`))) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`1p`))) +
  theme_legend +
  ggtitle("Chr 1p") +
  scale_color_manual(values = c("-1" = "blue",
                                "1" = "red",
                                "0" = "grey", 
                                "na" = "grey"))
ggsave("CNV_1p_leg.pdf", width = 10, height = 7)
ggsave("CNV_1p_noleg.pdf", width = 7, height = 7)


# Chr 14q
finaldf_CNV$`14q`[is.na(finaldf_CNV$`14q`)] <- "na"
FF <- subset(finaldf_CNV, finaldf_CNV$`14q` == 0) 
FF2 <- subset(finaldf_CNV, finaldf_CNV$`14q` == 1) 
FF3 <- subset(finaldf_CNV, finaldf_CNV$`14q` == -1) 
ggplot(
  finaldf_CNV %>% arrange(factor(`14q`)),
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(`14q`)
  )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`14q`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`14q`))) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`14q`))) +
  theme_nolegend +
  ggtitle("Chr 14q") +
  scale_color_manual(values = c("-1" = "blue",
                                "1" = "red",
                                "0" = "grey", 
                                "na" = "grey"))
ggsave("CNV_14q_leg.pdf", width = 10, height = 7)
ggsave("CNV_14q_noleg.pdf", width = 7, height = 7)


# Chr 22q
finaldf_CNV$`22q`[is.na(finaldf_CNV$`22q`)] <- "na"
FF <- subset(finaldf_CNV, finaldf_CNV$`22q` == 0) 
FF2 <- subset(finaldf_CNV, finaldf_CNV$`22q` == 1) 
FF3 <- subset(finaldf_CNV, finaldf_CNV$`22q` == -1) 
ggplot(
  finaldf_CNV %>% arrange(factor(`22q`)),
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(`22q`)
  )
) + 
  geom_point(size=1.5) +
  geom_point(data=FF, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`22q`))) +
  geom_point(data=FF2, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`22q`))) +
  geom_point(data=FF3, aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(`22q`))) +
  theme_nolegend +
  ggtitle("Chr 22q") +
  scale_color_manual(values = c("-1" = "blue",
                                "1" = "red",
                                "0" = "grey", 
                                "na" = "grey"))
ggsave("CNV_22q_leg.pdf", width = 10, height = 7)
ggsave("CNV_22q_noleg.pdf", width = 7, height = 7)


##
#Generate UMAPs for all chr arms
##

# Create a PDF file for output
pdf("SupplementalTable_CNAbyChrArm_1298_2.pdf")

# Initialize variables
plots_per_page <- 16  # 4x4 grid
total_plots <- ncol(CNVcalls2) - 1  # Excluding the first column

# Loop through plots
for (i in seq(1, total_plots, plots_per_page)) {
  end <- min(i + plots_per_page - 1, total_plots)
  
  # Create a list to hold plots for this page
  plot_list <- list()
  
  # Loop through plots on this page
  for (z in i:end) {
    ggdf <- cbind(finaldf2, cell = CNVcalls2[, z + 1])  # Skip first column
    title <- colnames(CNVcalls2)[z + 1]
    grey = subset(ggdf, ggdf$cell == "0")
    blue = subset(ggdf, ggdf$cell == "-1")
    red = subset(ggdf, ggdf$cell == "1")
    p <- ggplot(ggdf %>% 
                  arrange(factor(cell)),
                aes(
                  x = UMAP1_2D, 
                  y = UMAP2_2D,
                  colour = factor(cell)
                )) +
      geom_point(size = 0.8) +
      geom_point(data=grey,
                 aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(cell)), size = 0.8)+
      geom_point(data=red,
                 aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(cell)), size = 0.8)+
      geom_point(data=blue,
                 aes(x=UMAP1_2D, y=UMAP2_2D, color=factor(cell)), size = 0.8)+
      theme_blank() +
      theme(legend.position = "none",  # Remove individual legends
            axis.text = element_blank(),  # Remove axis text
            axis.title = element_blank()) +  # Remove axis titles
      ggtitle(title) +
      scale_color_manual(values = c("-1" = "blue",
                                    "1" = "red",
                                    "0" = "grey80"))
    
    plot_list[[z - i + 1]] <- p
  }
  # Combine plots for this page and arrange in a grid
  page_plots <- do.call(grid.arrange, c(plot_list, ncol = 4))
  # Create a combined legend for this page
  combined_legend <- get_legend(plot_list[[1]])  # Assuming all plots have the same legend
  # Arrange the plots and legend on the page
  page <- grid.arrange(arrangeGrob(page_plots, combined_legend,
                                   ncol = 2, widths = c(4, 1)))
  # Print the combined plot for this page
  print(page)
}
# Save and close the PDF file
dev.off()


####
#Kaplan Miere plots on tumors with and without copy number alterations CNA 
####

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
mycol <- c("Black", "Blue")

sample_info <-readxl::read_xlsx("~/DAM_data_clinical_patient_v40_20230808a.xlsx", skip=5)
sample_info <- sample_info[-1,]
log2tpm <- readRDS("~/1298samples_combatseq_log2tpm.rds")
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

# KM plot for samples with and without 10q loss in cluster A
CNVcalls <- read.csv("~/1298_CNVcalls.csv", check.names = FALSE)
chr10qLoss <- subset(CNVcalls, CNVcalls$`10q` == "-1")
timerec$chr10qLoss[timerec$coordinate_ID %in% chr10qLoss$coordinate_ID] <- "Yes"
timerec$chr10qLoss[!(timerec$coordinate_ID %in% chr10qLoss$coordinate_ID)] <- "No"
#select only tumors in cluster A
timerec_A <- subset(timerec, timerec$coordinate_ID %in% cluster_A$coordinate_ID)

fitAll <- surv_fit(Surv(time, status) ~ chr10qLoss, data = timerec_A)
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)",ggtheme = theme_KM , palette = mycol) 
ggsave("Chr10qLoss_outcome_leg.pdf", width = 7, height = 7)

# KM plot for samples with and without 6q loss 
CNVcalls <- read.csv("~/1298_CNVcalls.csv", check.names = FALSE)
chr6qLoss <- subset(CNVcalls, CNVcalls$`16q` == "-1")
timerec$chr6qLoss[timerec$coordinate_ID %in% chr6qLoss$coordinate_ID] <- "Yes"
timerec$chr6qLoss[!(timerec$coordinate_ID %in% chr6qLoss$coordinate_ID)] <- "No"
#select only tumors in cluster A
timerec_A <- subset(timerec, timerec$coordinate_ID %in% cluster_A$coordinate_ID)

fitAll <- surv_fit(Surv(time, status) ~ chr6qLoss, data = timerec_A)
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)",ggtheme = theme_KM , palette = mycol) 
ggsave("Chr6qLoss_outcome_leg.pdf", width = 7, height = 7)


# KM plot for samples with and without 22q loss 
CNVcalls <- read.csv("~/1298_CNVcalls.csv", check.names = FALSE)
chr22qLoss <- subset(CNVcalls, CNVcalls$`22q` == "-1")
timerec$chr22qLoss[timerec$coordinate_ID %in% chr22qLoss$coordinate_ID] <- "Yes"
timerec$chr22qLoss[!(timerec$coordinate_ID %in% chr22qLoss$coordinate_ID)] <- "No"
#select only tumors in cluster A
timerec_A <- subset(timerec, timerec$coordinate_ID %in% cluster_A$coordinate_ID)
fitAll <- surv_fit(Surv(time, status) ~ chr22qLoss, data = timerec)
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)",ggtheme = theme_KM , palette = mycol) 
ggsave("Chr22qLoss_outcome_leg.pdf", width = 7, height = 7)

# KM plot for samples with and without 14q loss 
CNVcalls <- read.csv("~/1298_CNVcalls.csv", check.names = FALSE)
chr14qLoss <- subset(CNVcalls, CNVcalls$`14q` == "-1")
timerec$chr14qLoss[timerec$coordinate_ID %in% chr14qLoss$coordinate_ID] <- "Yes"
timerec$chr14qLoss[!(timerec$coordinate_ID %in% chr14qLoss$coordinate_ID)] <- "No"
#select only tumors in cluster A
timerec_A <- subset(timerec, timerec$coordinate_ID %in% cluster_A$coordinate_ID)

fitAll <- surv_fit(Surv(time, status) ~ chr14qLoss, data = timerec_A)
ggsurvplot(fitAll,
           conf.int = FALSE,
           data = timerec,
           pval = TRUE,
           pval.method = FALSE,
           risk.table = FALSE, 
           combine = TRUE, 
           ylab = "Recurrence-free Rate", 
           xlab = "Time (months)",ggtheme = theme_KM , palette = mycol) 
ggsave("Chr14qLoss_outcome_leg.pdf", width = 7, height = 7)


