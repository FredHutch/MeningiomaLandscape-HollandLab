#title: DBSCAN, k-means and hierarchical clustering 
#authors: Nayanga Thirimanne

library(dbscan)
library(ggplot2)
library(dplyr)

########
#DBSCAN on umap 3D coordinates
########

#load UMAP 3D coordinates
umapcord <- readRDS("~/UMAP_versions/1298samples_umap3D_041923.rds")

#run dbscan
cl <- dbscan(umapcord[,2:4], eps=0.5, minPts = 10)
clust <- cl$cluster

# Merge cluster IDs with UMAP 2D coordinates (for plotting purposes)
finaldf <- readRDS("~/1298samples_umap2D_041923.rds")
finaldf_g <- finaldf %>%
  mutate(cluster = clust) 

ggplot(
  finaldf_g,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(cluster)
  )
) + 
  geom_point(size=2) +
  scale_color_manual(values = c("1" = "green", 
                                "2" = "blue", 
                                "3" = "red", 
                                "4" = "black",
                                "5" = "brown",
                                "6" = "gold", 
                                "7" = "cyan", 
                                "8" = "purple", 
                                "9" = "pink"
                                ))
table(clust)


#######
#k-means clustering
#######

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(FactoMineR) # clustering algorithms & visualization
install.packages("factoextra")
library(factoextra)

df21 <- readRDS("~/1298samples_combatseq_log2tpm.rds")
df <- t(scale(df21))

## determining optimum number of clusters
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:30
# extract wss for 2-30 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

k2 <- kmeans(df, centers = 14, nstart = 25)
str(k2)
k2_cluster <- as.data.frame(k2$cluster) %>% rownames_to_column(., "coordinate_ID")
finaldf_50 <- left_join(finaldf, k2_cluster, by="coordinate_ID")

ggplot(
  finaldf_50,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(k2$cluster)
  )
) + 
  geom_point(size=2) +
  theme_nolegend +
  scale_color_manual(values = c("1" = "green", 
                                "2" = "blue", 
                                "3" = "red", 
                                "4" = "black",
                                "5" = "brown",
                                "6" = "gold", 
                                "7" = "cyan", 
                                "8" = "purple", 
                                "9" = "pink",
                                "10" = "grey", 
                                "11" = "yellow4", 
                                "12" = "orange", 
                                "13" = "goldenrod4",
                                "14" = "forestgreen"))

ggsave("kmeans_leg.pdf", width = 7, height = 7)
ggsave("kmeans_noleg.pdf", width = 7, height = 7)



###########################
#hierarchical clustering 
###########################

# Dissimilarity matrix
d <- dist(df, method = "euclidean")
# Ward's method
hc5 <- hclust(d, method = "ward.D" )
# Cut tree
sub_grp <- cutree(hc5, k = 10)
# Number of members in each cluster
table(sub_grp)
# Plot the obtained dendrogram
plot(hc5, cex = 0.6, hang = -1)

subgrp <- as.data.frame(sub_grp)
subgrp <- rownames_to_column(subgrp, var = "coordinate_ID")
finaldf_e <- left_join(finaldf, subgrp, by="coordinate_ID")

ggplot(
  finaldf_e,
  aes(
    x = UMAP1_2D, 
    y = UMAP2_2D,
    color = factor(sub_grp)
  )
) + 
  geom_point(size=2) +
  theme_nolegend +
  scale_color_manual(values = c("1" = "green", 
                                "2" = "blue", 
                                "3" = "red", 
                                "4" = "black",
                                "5" = "brown",
                                "6" = "gold", 
                                "7" = "cyan", 
                                "8" = "purple", 
                                "9" = "pink",
                                "10" = "grey"))

ggsave("hclust_leg.pdf", width = 7, height = 7)
ggsave("hclust_noleg.pdf", width = 7, height = 7)

