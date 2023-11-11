# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)

# ===== Initialization =====
metadata <- readRDS("/data/metadata.RDS")
umap_embeddings <- readRDS("/data/umap_embeddings_Tcells.rda")

# ===== Figure 4a =====

metadata <- metadata[metadata$Compartment %in% c("CD8 T Cells","CD4 T Cells","NK T Cells"),]
setDT(metadata)
umap_embeddings$Cell.id <- rownames(umap_embeddings)
umap_embeddings$Alias <- NULL
setDT(umap_embeddings)

umap_embeddings <- merge(umap_embeddings,metadata, by="Cell.id")

sum_cdr3 <- umap_embeddings[,.N, by = .(orig.ident, cdr3s_aa)]

no_tcells <- vector()
Tcells <- list()

for (l in unique(sum_cdr3$orig.ident)){
  no_tcells <- c(no_tcells, sum(!is.na(umap_embeddings[orig.ident == l & Compartment %in% c("CD8 T Cells","CD4 T Cells","NK T Cells")]$cdr3s_aa)))
}


Tcells <- setnames(data.table(unique(sum_cdr3$orig.ident),no_tcells), c("orig.ident", "no_tcells"))

plot_df <- merge(umap_embeddings, sum_cdr3, by = c("orig.ident","cdr3s_aa"), all=T)
plot_df <- unique(plot_df)
plot_df[is.na(`cdr3s_aa`)]$`N` <- 0
plot_df$cdr3s_aa <- NULL

plot_df <- plot_df[,noCdr3s_aa := sum(N), by = orig.ident]
plot_df <- merge(plot_df, Tcells, by = "orig.ident")
plot_df$`Clonesize_proportion` <- (plot_df$N/plot_df$no_tcells)*100

plot_df$Type <- "T Cells"
plot_df[plot_df$N == 1]$Type <- "TCR-clonotype freq = 1"
plot_df[plot_df$N > 1]$Type <- "TCR-clonotype freq > 1"
plot_df[plot_df$N > 4]$Type <- "TCR-clonotype freq > 4"


plot_df$Type <- factor(plot_df$Type, levels = c("T Cells","TCR-clonotype freq = 1", "TCR-clonotype freq > 1", "TCR-clonotype freq > 4"))

setDF(plot_df)
rownames(plot_df) <- plot_df$sample
plot_df$TimePoint_Final <- factor(plot_df$TimePoint_Final , levels = c("Leukapheresis", 
                                                                       "4 weeks after CAR T","6 months after CAR T",
                                                                       "12 months after CAR T"))

plot_df$Final_Alias <- plot_df$Alias
plot_df[plot_df$Compartment != "T Cells",]$Final_Alias <- plot_df[plot_df$Compartment != "T Cells",]$Compartment

plot_df <- plot_df[,c("TimePoint_Final","UMAP_1", "UMAP_2","N","Clonesize_proportion","Alias","Type","Patient_ID")] 
setnames(plot_df,"N","Clonosize")
cells <-sum(no_tcells)
plot_df1 <- plot_df[,c("UMAP_1", "UMAP_2","Alias")] 
plot_df2 <- plot_df1 %>% group_by(Alias) %>% summarize_all(mean)

p <-  ggplot(plot_df,  aes_string(x = "UMAP_1", y = "UMAP_2", color="Type")) +
  geom_point(data = plot_df, pch = 21, size = 0.005,alpha = 0.2) +
  ggtitle(paste0("a/b TCRs (n=",cells,")")) + 
  scale_color_manual(values = c("grey","blue","#00bebe","yellow"), 
                     labels = paste0(levels(factor(plot_df$Type)))) +
  ggrepel::geom_text_repel(data = plot_df2,aes(label = Alias), size=2.5, max.iter = 60, color="black") +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=10, face='bold', hjust = 0.08),
        axis.title=element_text(size=8,color="black",face="bold"),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        legend.text = element_text(colour="black", size = 8),
        legend.key = element_rect(fill = "transparent", size = 12,colour = "transparent"),
        legend.key.size = unit(0.2, "cm"),  legend.position = c(0.85,1),
        legend.key.width = unit(1,"cm"), 
        axis.line = element_line(size = 0.5, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=1.5, alpha = 2)))

print(p)

