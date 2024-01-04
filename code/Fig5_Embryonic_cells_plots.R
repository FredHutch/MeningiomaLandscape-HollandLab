#title: Generate plots to visualize mouse embryonic cell types enriched in each meningioma cluster
#authors: Nayanga Thirimanne

library(dplyr)
library(ggplot2)

#load gene module score for each cluster
dat = readRDS('~/Analysis33_mouse_embryonic/data.rds')
gene_module_list = colnames(dat)[1:10]

major_trajectory_color_plate = c("Neuroectoderm_and_glia"             = "#f96100",
                                 "Intermediate_neuronal_progenitors"  = "#2e0ab7",
                                 "Eye_and_other"             = "#00d450",
                                 "Ependymal_cells"           = "#b75bff",
                                 "CNS_neurons"               = "#e5c000",
                                 "Mesoderm"                  = "#bb46c5",
                                 "Definitive_erythroid"      = "#dc453e",
                                 "Epithelium"                = "#af9fb6",
                                 "Endothelium"               = "#00a34e",
                                 "Muscle_cells"              = "#ffa1f5",
                                 "Hepatocytes"               = "#185700",
                                 "White_blood_cells"         = "#7ca0ff",
                                 "Neural_crest_PNS_glia"     = "#fff167",
                                 "Adipocytes"                = "#7f3e39",
                                 "Primitive_erythroid"       = "#ffa9a1",
                                 "Neural_crest_PNS_neurons"  = "#b5ce92",
                                 "T_cells"                   = "#ff9d47",
                                 "Lung_and_airway"           = "#02b0d1",
                                 "Intestine"                 = "#ff007a",
                                 "B_cells"                   = "#01b7a6",
                                 "Olfactory_sensory_neurons" = "#e6230b",
                                 "Cardiomyocytes"            = "#643e8c",
                                 "Oligodendrocytes"          = "#916e00",
                                 "Mast_cells"                = "#005361",
                                 "Megakaryocytes"            = "#3f283d",
                                 "Testis_and_adrenal"        = "#585d3b")


### Plot embryonic cell types enriched in each meningioma cluster (median, horizontal plot, without legend) ###

for(i in 1:10){
  name = gene_module_list[i]
  print(name)
  
  df = dat[,c("major_cell_cluster", "celltype")]
  df$score = as.vector(dat[[name]])
  
  df_x = df %>% group_by(major_cell_cluster, celltype) %>% summarise(median_score = median(score))
  df_x = df_x %>%
    group_by(major_cell_cluster) %>%
    arrange(median_score, .by_group=TRUE)
  
  df$celltype = factor(df$celltype, levels = as.vector(df_x$celltype))
  
  try(ggplot(df, aes(celltype, score, fill = major_cell_cluster)) + geom_boxplot(outlier.shape = NA) + 
        labs(x="celltype", y="gene_module_score", title=name) +
        theme_classic(base_size = 12) +
        scale_fill_manual(values=major_trajectory_color_plate) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_blank(), 
              axis.text.y = element_text(color="black", size = 12)) +
        ggsave(paste0(name, ".pdf"),
               dpi = 300,
               height  = 6, 
               width = 30))
}

### Plot major cell clusters enriched in each meningioma cluster ###

dat2 <- dat[, c(1:7,11)]
datmean <- dat2 %>% group_by(major_cell_cluster) %>%
  summarise(across(everything(), mean),
            .groups = 'drop') %>% 
  as.data.frame()

datmeanlong <- melt(datmean)

legend_pt_size = 4
axis_text_size = 25
axis_title_size = 25
legend_text_size = 15
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_legend_em <- theme_classic()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank())

ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  labs(x=" ", y="gene module score") +
  theme_legend_em +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black", size = 15 , angle = 90), axis.text.y = element_text(color="black" , size = 15), 
        axis.title = element_text(size = 15, face = "bold"))

ggsave("Clusters_embryocelltype_mean_leg.pdf", width = 15, height = 7)

ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  labs(x=" ", y="gene module score") +
  theme_legend_em +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black", size = 15 , angle = 90), axis.text.y = element_text(color="black" , size = 15), 
        axis.title = element_text(size = 15, face = "bold")) +
  theme(legend.position = "none")
ggsave("Clusters_embryocelltype_mean_noleg.pdf", width = 10, height = 7)

