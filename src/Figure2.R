# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(pheatmap)
library(dplyr)
library(scater)
library(Seurat)
library(plyr)

# ===== Functions =====

source("/src/pheatmap_updated.R")

data_summary <- function(data, varname="Clonosize", 
                         groupnames=c("Condition","Alias")){
  require(plyr)
  se = function(x, na.rm=T){sd(x, na.rm=na.rm)/sqrt(length(x))}
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


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
metadata <- readRDS("/data/metadata.RDS")
cytotoxic_activity_scores <- fread("/data/cytotoxic_activity_scores.RDS")
umap_embeddings <- readRDS("/data/umap_embeddings_Tcells.rda")
counts <- readRDS("/data/normalized_count_data.RDS")
markers_tcells <- fread("/data/markers_Tcells.txt", sep="\t")


tissue_so <- CreateSeuratObject(counts = counts, project = "ZUMA1", assay = "RNA",
                                min.cells = 0, min.features = 0, names.field = 1,
                                names.delim = "-", meta.data = metadata)


# ===== Figure 1a =====

###### Remove samples with Unknown response and cells annotated as doublets. Retain only CD4/CD8 T Cells.

setDT(metadata)
metadatametadata <- metadata[!metadata$orig.ident %in% c("R1","R2","T8"),]
metadata <- metadata[Response != "Unknown"]
metadata <- metadata[!Alias %in% c("DBL1","DBL2","DBL3","DBL4","DBL5","DBL6","DBL7","DBL8","DBL9")]
metadata <- metadata[!(Patient_ID == "P10" & TimePoint_Final == "12 months after CAR T")]
metadata$Label <- paste0(metadata$TimePoint_Final,"_",metadata$orig.ident,"@", metadata$Response)

metadata$Compartment.new <- metadata$Compartment
metadata[regexpr("CD8.*",metadata$Alias)==TRUE]$Compartment.new <- "CD8 T cells"
metadata[regexpr("CD4.*|T regs",metadata$Alias)==TRUE]$Compartment.new <- "CD4 T cells"

metadata <- metadata[Compartment.new %in% c("CD8 T cells","CD4 T cells")]
metadata$Condition <- metadata$Response

sampleFreq <- setnames(plyr::count(metadata$Label), c("Sample", "All_Cells"))

###### Calculate CD4 and CD8 proportions

total_proprotions_across_groups.1 <- setnames(setDT(plyr::count(metadata[,c("Label","TimePoint_Final","Condition","Compartment.new"), with=F])),
                                              c("Sample","TimePoint_Final","Condition","Compartment.new","Counts"))

total_proprotions_across_groups <- dcast(total_proprotions_across_groups.1, Sample+TimePoint_Final+Condition~Compartment.new)

total_proprotions_across_groups[is.na(total_proprotions_across_groups)] <- 0

total_proprotions_across_groups.1 <- melt(total_proprotions_across_groups)
setnames(total_proprotions_across_groups.1, c("variable","value"), c("Compartment.new","Counts"))

total_proprotions_across_groups.1 <- merge(total_proprotions_across_groups.1, sampleFreq, by="Sample")

total_proprotions_across_groups.1$Cellular_Proportion = (total_proprotions_across_groups.1$Counts/total_proprotions_across_groups.1$All_Cells)*100


total_proprotions_across_groups <- setDT(total_proprotions_across_groups)
total_proprotions_across_groups.1 <- setDT(total_proprotions_across_groups.1)


total_proprotions_across_groups.1 <- setnames(setDT(plyr::count(metadata[,c("Label","TimePoint_Final","Condition","Compartment.new"), with=F])),
                                              c("Sample","TimePoint_Final","Condition","Compartment.new","Counts"))

total_proprotions_across_groups <- dcast(total_proprotions_across_groups.1, Sample+TimePoint_Final+Condition~Compartment.new)

total_proprotions_across_groups[is.na(total_proprotions_across_groups)] <- 0

total_proprotions_across_groups$CD8.CD4.ratio <- total_proprotions_across_groups$`CD8 T cells`/total_proprotions_across_groups$`CD4 T cells`
total_proprotions_across_groups$`CD4 T cells` <- NULL
total_proprotions_across_groups$`CD8 T cells` <- NULL

total_proprotions_across_groups.1 <- melt(total_proprotions_across_groups)
setnames(total_proprotions_across_groups.1, c("variable","value"), c("Compartment.new","CD8.CD4.ratio"))

setDT(total_proprotions_across_groups.1)

###### Plot CD8/CD4 ratio

df2 <-total_proprotions_across_groups.1
df2$Condition <- factor(df2$Condition, levels = c("Long-term Responder","Relapsed","non-Responder"))
df2$Alias <- df2$Compartment.new

colors <- c("blue3","darkolivegreen3","brown3")

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4 weeks after CAR T", "6 months after CAR T", "12 months after CAR T"))
df2$Condition <- factor(df2$Condition, levels = c("Long-term Responder","Relapsed","non-Responder"))

y.value <- max(df2$CD8.CD4.ratio)

setDT(df2)


ratio.plot <- ggplot(data = df2, aes(x=`Condition`, y=`CD8.CD4.ratio`, fill=Condition)) + facet_wrap(~`TimePoint_Final`, ncol = 2) +  
  geom_boxplot(data = df2,alpha=0.6, aes(x=`Condition`, y=`CD8.CD4.ratio`), width=0.3, outlier.size = FALSE,position=position_dodge(width = 0.9))+ 
  ylab("CD8/CD4 T Cell ratio")+ xlab("")+ ggtitle("") + 
  ylim(0,(y.value+2)) +
  geom_jitter(data = df2, aes( x = Condition, y =CD8.CD4.ratio),position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=0.5,fill="grey") + 
  scale_fill_manual(values = colors) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=8,color="black",face="bold"), 
        axis.title=element_text(size=9,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =9),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position ='top', 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.x = element_blank(),
        strip.text = element_text(size=9)) 


print(ratio.plot)


# ===== Figure 1b =====

cytotoxic_activity_scores$Response <- factor(cytotoxic_activity_scores$Response, levels=c("Long-term Responder","Relapsed","non-Responder"))
cytotoxic_activity_scores$TimePoint_Final <- factor(cytotoxic_activity_scores$TimePoint_Final, levels=c("Leukapheresis","4 weeks after CAR T","6 months after CAR T", "12 months after CAR T"))

y.value <- max(cytotoxic_activity_scores$`cytotoxic_activity_UCell`)
dodge <- position_dodge(width = 0.8)
colors <- c("blue3","darkolivegreen3","brown3")

signature <- ggplot(data = cytotoxic_activity_scores, aes(x=`TimePoint_Final`, y=`cytotoxic_activity_UCell`, fill=Response)) +  
  ylab("Cytotoxic Signature")+ xlab("")+ 
  geom_violin(aes(fill = Response), width=0.8, position=dodge,alpha=0.6) +
  geom_boxplot(aes(group=interaction(Response,TimePoint_Final)), 
               width=0.2, fill="white", position=dodge,
               outlier.shape=NA)+scale_fill_manual(values = colors) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=12,color="black",face="bold"), 
        axis.title=element_text(size=10,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =10),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = 'top', 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=10, face="bold"),
        strip.text = element_text(size=10))

print(signature)

# ===== Figure 1c =====

###### UMAP with T cell populations

umap_embeddings$Alias <- as.factor(umap_embeddings$Alias )

# Make color scale
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

plot_df <- umap_embeddings %>% group_by(Alias) %>% summarize_all(mean)

umapTcells <- ggplot(umap_embeddings, aes(x = UMAP_1, y = UMAP_2)) + ggtitle(paste0("T cells, n=", nrow(umap_embeddings))) + 
  geom_point(data = umap_embeddings, pch = 21, color = 'black', size = 0.005,alpha = 0.4) +
  geom_point(aes(color = Alias), alpha = 0.4, size = 0.005, show.legend = T) +
  ggrepel::geom_text_repel(data = plot_df,aes(label = Alias), size=4, max.iter = 80) +
  scale_color_manual(values = cols, 
                     labels = paste0(levels(factor(umap_embeddings$Alias)))) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust = 0.08),
        axis.title=element_text(size=8,color="black",face="bold"),
        legend.title = element_blank(),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        legend.text = element_text(colour="black", size = 6),
        legend.key = element_rect(fill = "transparent", size = 2,colour = "transparent"),
        legend.key.size = unit(1, "cm"), legend.position = "none", 
        axis.line = element_line(size = 0.5, colour = "black"))


print(umapTcells)


###### Heatmap with CD4 markers

subset_so = tissue_so[,regexpr("CD4.*|T regs",tissue_so$Alias)==TRUE]

markers_cd4 <- markers_tcells[markers_tcells$Category == "CD4",]
markers_cd4 <- unique(markers_cd4)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

SubClusterAvg = sumCountsAcrossCells(
  x = subset_so@assays$RNA@data,
  ids = subset_so$Alias,
  average = T)

SubClusterAvg <- SubClusterAvg[,colnames(SubClusterAvg) %in% subset_so$Alias]

dt <- as.data.frame(SubClusterAvg[rownames(SubClusterAvg) %in% as.character(markers_cd4$Marker),])
dt$Marker <- rownames(dt)
dt <- merge(dt, markers_cd4, by = "Marker")
dt <- setDF(dt)
rownames(dt) <- dt$Marker 
dt$Marker <- NULL
dt$Gene <- NULL
dt$Category <- NULL

dt <- na.omit(dt)
dt <- scale_mat(dt, "row")

myBreaks <- c(seq(min(dt), 0, length.out=3), 
              seq(1, max(dt), length.out=6))

heatmapCD4<- pheatmap_updated(mat = dt,cluster_rows=TRUE,cluster_cols = TRUE,scale="column",
                    breaks =myBreaks,
                    color =inferno(10),angle_col = "45",fontsize = 7,
                    name="Log(Expr)",heatmap_legend_param = list(title_gp=gpar(fontsize=12,fontface="bold"),labels_gp = gpar(fontsize = 12)))


print(heatmapCD4)


###### Heatmap with CD8/NK T markers

subset_so = tissue_so[,regexpr("CD8.*|NK T.*",tissue_so$Alias)==TRUE]

markers_cd8 <- markers_tcells[markers_tcells$Category == "CD8",]
markers_cd8 <- unique(markers_cd8)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

SubClusterAvg = sumCountsAcrossCells(
  x = subset_so@assays$RNA@data,
  ids = subset_so$Alias,
  average = T)

SubClusterAvg <- SubClusterAvg[,colnames(SubClusterAvg) %in% subset_so$Alias]

dt <- as.data.frame(SubClusterAvg[rownames(SubClusterAvg) %in% as.character(markers_cd8$Marker),])
dt$Marker <- rownames(dt)
dt <- merge(dt, markers_cd8, by = "Marker")
dt <- setDF(dt)
rownames(dt) <- dt$Marker 
dt$Marker <- NULL
dt$Gene <- NULL
dt$Category <- NULL

dt <- na.omit(dt)
dt <- scale_mat(dt, "row")

myBreaks <- c(seq(min(dt), 0, length.out=3), 
              seq(1, max(dt), length.out=6))

heatmapCD8<- pheatmap_updated(mat = dt,cluster_rows=TRUE,cluster_cols = TRUE,scale="column",
                              breaks =myBreaks,fontsize = 7,
                              color =inferno(10),angle_col = "45",
                              name="Log(Expr)",heatmap_legend_param = list(title_gp=gpar(fontsize=12,fontface="bold"),labels_gp = gpar(fontsize = 12)))


print(heatmapCD8)
