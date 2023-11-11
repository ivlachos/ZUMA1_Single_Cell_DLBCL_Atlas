# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(cowplot)
library(pheatmap)
library(plyr)
library(grid)
library(gridExtra)
library(ggplotify)

# ===== Functions =====

plot_no_of_clonotypes <- function(sum_cdr3, name){
  
  clons2 <- dcast(sum_cdr3, cdr3s_aa ~ Alias, sum, value.var = "N")
  
  total_df <- data.frame()
  for (i in seq(2,length(colnames(clons2[2:ncol(clons2)])))){
    for (j in seq(2,length(colnames(clons2[2:ncol(clons2)])))){
      x <- clons2[,i, with=F]
      y <- clons2[,j, with=F]
      int_df <- data.frame(x, y)
      setnames(int_df, c("x","y"))
      int_df <- setDT(int_df)
      sumTotal <- nrow(int_df[x>0 & y>0])
      total_df[i-1,j-1] <- sumTotal
    }
  }
  colnames(total_df) <- colnames(clons2[,2:ncol(clons2)])
  rownames(total_df) <- colnames(clons2[,2:ncol(clons2)])
  
  total_df <- log(total_df+1, base=2)
  
  p1<- pheatmap(mat = total_df,cluster_rows=TRUE,cluster_cols=TRUE, main = paste0(name),scale="none",
                fontsize = 10,fontsize_row =10, fontsize_col = 10,width = 10,height = 10)
  
  return(p1)
}

# ===== Initialization =====
metadata <- readRDS("/data/metadata.RDS")

# ===== Figure 4c =====

metadata <- setDT(metadata)
metadata <- metadata[Compartment %in% c("CD8 T Cells","CD4 T Cells","NK T Cells") & !is.na(cdr3s_aa)]
metadata <- metadata[Response == "Long-term Responder"]

clono_seurat_subset <- unique(metadata[,c("Patient_ID","cdr3s_aa","orig.ident","TimePoint_Final","Alias"), with=F])

sum_cdr3 <- clono_seurat_subset[!is.na(clono_seurat_subset$cdr3s_aa),]
sum_cdr3 <- sum_cdr3[,.N, by = .(TimePoint_Final, Alias,orig.ident, cdr3s_aa)]

clono_seurat_subset_noTcells <- setnames(setDT(plyr::count(clono_seurat_subset[!is.na(clono_seurat_subset$cdr3s_aa)][,c("orig.ident","TimePoint_Final","Patient_ID"), with=F])),
                                         c("orig.ident", "TimePoint_Final", "noTcells","Patient_ID"))

sum_cdr3 <- merge(sum_cdr3, clono_seurat_subset_noTcells, by = c("orig.ident","TimePoint_Final"))


p1 <- plot_no_of_clonotypes( sum_cdr3[TimePoint_Final == "Leukapheresis"], "Leukapheresis")
p2 <- plot_no_of_clonotypes( sum_cdr3[TimePoint_Final == "4 weeks after CAR T"], "4 weeks after CAR T")
p3 <- plot_no_of_clonotypes( sum_cdr3[TimePoint_Final == "6 months after CAR T"], "6 months after CAR T")
p4 <- plot_no_of_clonotypes( sum_cdr3[TimePoint_Final == "12 months after CAR T"], "12 months after CAR T")


plot_grid(plotlist = list(as.grob(p1),as.grob(p2),as.grob(p3),as.grob(p4)),  ncol = 2, nrow=2)
