# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(scater)

# ===== Functions =====

source("/src/pheatmap_updated.R")

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

# ===== Initialization =====

counts <- readRDS("/data/normalized_count_data.RDS")
metadata <- readRDS("/data/metadata.RDS")
umap_embeddings <- readRDS("/data/umap_embeddings.RDS")
major_markers <- fread("/data/major_markers.txt", sep="\t")
major_markers_heatmap <- fread("/data/major_markers_heatmap.txt", sep="\t")


tissue_so <- CreateSeuratObject(counts = counts, project = "ZUMA1", assay = "RNA",
                                min.cells = 0, min.features = 0, names.field = 1,
                                names.delim = "-", meta.data = metadata)

# ===== Figure 1c =====

umap_embeddings$Alias <- as.factor(umap_embeddings$Alias)

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

plot_df <- umap_embeddings %>% group_by(Compartment) %>% summarize_all(mean)


figureA <- ggplot(umap_embeddings, aes(x = UMAP1, y = UMAP2)) + ggtitle("n=405,775") + 
        geom_point(data = umap_embeddings, pch = 21, color = 'black', size = 0.005,alpha = 0.4) +
        geom_point(aes(color = Alias), alpha = 0.4, size = 0.005, show.legend = T) +
        ggrepel::geom_text_repel(data = plot_df,aes(label = Compartment), size=3.5, max.iter = 70) +
        scale_color_manual(values = cols,labels = paste0(levels(factor(umap_embeddings$Alias)))) +
        theme(panel.border = element_blank(), panel.background = element_blank(),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              plot.title = element_text(size=10, face='bold', hjust = 0.08),
              axis.title=element_text(size=8,color="black",face="bold"),
              legend.title = element_blank(),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
              legend.text = element_text(colour="black", size = 6),
              legend.key = element_rect(fill = "transparent", size = 2,colour = "transparent"),
              legend.key.size = unit(1, "cm"), legend.position = "none", 
              axis.line = element_line(size = 0.5, colour = "black"))


print(figureA)


# ===== Figure 1d =====

cellEmbed <- umap_embeddings

plot <- list();
k <- 1
for(j in major_markers$Gene_rnaSeq){
  
  assay_data <- GetAssayData(object = tissue_so, assay= "RNA")[j,]
  df = data.frame(
    x=cellEmbed[, 1], 
    y=cellEmbed[, 2])
  
  df$`Log2(Expr)`=assay_data
  
  df$alpha <- 1
  df[df$`Log2(Expr)`!=0,"alpha"] <- 0.5
  df$size <- 0.005
  df[df$`Log2(Expr)`!=0,"size"] <- 2.5
  data<-df[order(df$`Log2(Expr)`, decreasing=FALSE),]
  
  
  name <- paste0(major_markers[Gene_rnaSeq == j]$Gene_rnaSeq, ", ",   major_markers[Gene_rnaSeq == j]$Type)  
  plot[[k]] <- ggplot(data,aes(x=x, y=y, colour=`Log2(Expr)`),mar=c(0,0,3,0)) + 
    ggtitle(name) +
    geom_point(data = data, pch = 21, size = 0.00005,alpha = 0.3) + 
    ylab("UMAP2") + xlab("UMAP1") + 
    scale_color_gradientn(colours = colorRampPalette(c("darkblue","yellow"))(3)) + #scale_size(range=c(0,0.7))+
    
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(size=6,color="black",face="bold"), 
          axis.title=element_blank(),
          legend.title = element_text(size=6,color="black",face="bold"),
          legend.text = element_text(colour="black", size = 5),
          legend.key = element_rect(fill = "transparent", size = 5,colour = "transparent"),
          legend.key.size = unit(0.2, "cm"), legend.direction="horizontal", 
          legend.key.width = unit(0.2,"cm"), legend.position = c(0.6,0.9999),
          axis.line = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")) 
  
  
  min_y <- min(cellEmbed[,2])
  min_x <- min(cellEmbed[,1])
  max_y <- max(cellEmbed[,2])
  max_x <- max(cellEmbed[,1])
  
  statement1 <- 'atop(bold("UMAP1"))'
  statement2 <- 'atop(bold("UMAP2"))'
  
  
  plot[[k]] <- plot[[k]] + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y)) + 
    geom_segment(aes(y = min_y,
                     yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                     x = min_x - 1,
                     xend = min_x - 1),
                 arrow = arrow(length = unit(0.1, "cm")), color="black")+
    geom_segment(aes(y = min_y,
                     yend = min_y,
                     x = min_x - 1,
                     xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1), color="black",
                 arrow = arrow(length = unit(0.1, "cm"))) +
    annotate(geom="text", x=min_x + 5.5, y=min_y - 2,
             color="black", label = statement1,size= 2, parse = TRUE) + 
    annotate(geom="text", x=min_x + 0.5, y=min_y + 7,color="black", label = statement2,size= 2, parse = TRUE,angle = 90)
  
  
  k <- k+1
  
}

plot_grid(plotlist = plot,ncol = 3, nrow=6)


# ===== Figure 1e =====


SubClusterAvg = sumCountsAcrossCells(
  x = tissue_so@assays$RNA@data,
  ids = tissue_so$Major_Alias,
  average = T)

Major_Alias <- unique(tissue_so$Major_Alias)
Major_Alias <- Major_Alias[Major_Alias != "Doublets"]

SubClusterAvg <- SubClusterAvg[,colnames(SubClusterAvg) %in% Major_Alias]

dt <- as.data.frame(SubClusterAvg[rownames(SubClusterAvg) %in% as.character(major_markers_heatmap$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, major_markers_heatmap, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene_rnaSeq 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Type` <- NULL

dt <- na.omit(dt)


dt <- scale_mat(dt, "row")


figureC = draw(Heatmap(as.matrix(dt),
                   clustering_distance_columns = "euclidean",  # Set Euclidean distance for column clustering
                   clustering_distance_rows = "euclidean",  # Set Euclidean distance for row clustering
                   clustering_method_columns = "complete",  # Set linkage method, e.g., "complete"
                   clustering_method_rows = "complete",
                   column_names_rot = 50,  
                   cluster_columns  = TRUE,  
                   cluster_rows = TRUE,  
                   row_names_gp = gpar(fontsize = 11),name="Log(Expr)",
                   column_names_gp = gpar(fontsize = 11), 
                   show_column_names = TRUE,
                   heatmap_height = unit(0.6, "cm")*nrow(dt),
                   heatmap_width = unit(0.65, "cm")*ncol(dt),
                   
                   heatmap_legend_param = list(#title_gp = gpar(col = "black", 
                     title_gp = gpar(fontsize = 11, fontface = "bold"), 
                     fontsize = 25,face="bold",
                     legend_direction = "horizontal", 
                     legend_width= unit(6, "cm"),labels_gp = gpar(fontsize = 11))),
           
           
           heatmap_legend_side = "top")


print(figureC)


