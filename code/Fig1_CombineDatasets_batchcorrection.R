#title: Combine meningioma datasets, perform batch correction and normalize
#authors: Nayanga Thirimanne and Sonali Arora

library(sva)
library(dplyr)
library(rtracklayer)
library(edgeR)

##### load rawcounts from each dataset  #####
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

#####  batch correction #####
DF2 <- cbind(df1, df2, df3, df4, df5, df6, 
             df7, df8, df9, df10, df11, df12, df13, df14, df15)
DF2 <- as.data.frame(DF2)
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
           rep(15, ncol(df15))
)

adjusted_raw <- ComBat_seq(as.matrix(DF2), batch=batch, group=NULL)
adjusted_raw_2 <- as.data.frame(adjusted_raw)

##### select protein coding genes and normalize read counts  #####
# load gencode gtf file
gtf <- rtracklayer::import("~/Reference_Genomes/hg38_release39_053122/gencode.v39.primary_assembly.annotation.gtf")
hg38gtf <- as.data.frame(gtf)
#extract only important columns and protein coding genes
hg38gtfsub1 <- hg38gtf[ ,c(1:4,7,10:12)]
hg38gtfsub2 <- filter(hg38gtfsub1, hg38gtfsub1$gene_type == "protein_coding" & hg38gtfsub1$type == "gene")
#remove gene_id version number
hg38gtfsub2$gene_id <- gsub("\\.[0-9]*$", "", hg38gtfsub2$gene_id)
# select protein coding genes
AB <- left_join(adjusted_raw_2 %>%
                  mutate(gene_name = rownames(adjusted_raw_2)), hg38gtfsub2, by="gene_name", multiple = "all")
AB1 <- dplyr::filter(AB, AB$gene_type == "protein_coding") %>%
  arrange(gene_name) %>%
  dplyr::filter(duplicated(gene_name) == FALSE)
rownames(AB1) <- AB1$gene_name
drops <- c("gene_id", "seqnames", "start", "end", "strand", "gene_type", "type", "gene_name", "width")
AB2 <- AB1[ , !(names(AB1) %in% drops)]
# normalize
gene_width = AB1$width
gene_name = rownames(adjusted_raw_2)
y <- DGEList(counts = adjusted_raw_2)
rpkm = rpkm(y, gene.length = gene_width)
rownames(rpkm) = gene_name
tpm = apply(rpkm, 2, function(x){
  (x/sum(x))*10^6
})
log2_tpm_adjusted =log2(tpm+1) #batch corrected, normalized read counts

