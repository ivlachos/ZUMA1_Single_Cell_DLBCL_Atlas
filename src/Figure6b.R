# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(ggplot2)
library(viridis)
library(plyr)
library(grid) 

# ===== Initialization =====

metadata <- readRDS("/data/metadata.RDS")
metadata <- setDT(metadata)
metadata <- metadata[!Alias %in% c("DBL1","DBL2","DBL3","DBL4","DBL5","DBL6","DBL7","DBL8")]
metadata <- metadata[Major_Alias %in% c("CD8 T Cells","CD4 T Cells","NK T Cells") & !is.na(cdr3s_aa)]

# ===== Figure 6b =====

tcrs.4weeks <- fread("/data/TCRs_Leukapheresis_4 weeks_CytotoxicCombined.txt")

tcrs.total <- tcrs.4weeks[Response == "Long-term Responder"]
tcrs.total <- tcrs.4weeks[Alias == "T cells cytotoxic"]
tcrs.total <- tcrs.total[Type == "Expanded"]

clono_seurat_subset = metadata
clono_seurat_subset[Alias %in% c("CD8 Tem/Highly cytotoxic","CD8 Tem/Proinflammatory","CD8 Tem/TGFb1+CXCR4+","CD4 CTL")]$Alias <- "T cells cytotoxic"
clono_seurat_subset$TimePoint_Final <- factor(clono_seurat_subset$TimePoint_Final , levels = c("Leukapheresis", 
                                                                                               "4wk post-CAR T","6mo post-CAR T",
                                                                                               "12mo post-CAR T"))

clono_seurat_subset$id <- rownames(clono_seurat_subset)
clono_seurat_subset <- unique(clono_seurat_subset[,c("id","orig.ident","cdr3s_aa","orig.ident","TimePoint_Final")])
clono_seurat_subset <- setDT(clono_seurat_subset)

sum_cdr3 <- clono_seurat_subset[,.N, by = .(orig.ident,TimePoint_Final, cdr3s_aa)]
sum_cdr3 <- sum_cdr3[!is.na(sum_cdr3$cdr3s_aa),]

clono_seurat_subset_noTcells <- setnames(setDT(plyr::count(clono_seurat_subset[!is.na(clono_seurat_subset$cdr3s_aa)][,c("orig.ident","TimePoint_Final"), with=F])),
                                               c("orig.ident","TimePoint_Final", "noTcells"))


sum_cdr3 <- sum_cdr3[,noCdr3s_aa := sum(N), by = orig.ident]

sum_cdr3 <- merge(sum_cdr3, clono_seurat_subset_noTcells[,c("orig.ident","noTcells"), with=F], by = "orig.ident")

sum_cdr3$`Clonesize_proportion` <- (sum_cdr3$N/sum_cdr3$noTcells)*100


setupinfo <- unique(metadata[,c("orig.ident","Response"), with=F])

sum_cdr3 <- merge(sum_cdr3, setupinfo, by = "orig.ident")

top.clonotypes.responders <- sum_cdr3[order(-Clonesize_proportion),][Response=="LtR"][cdr3s_aa %in% tcrs.total[Response == "Long-term Responder"]$cdr3s_aa]


sum_cdr3.responders <- sum_cdr3[Response=="LtR"][cdr3s_aa %in% top.clonotypes.responders$cdr3s_aa]


sum_total <- sum_cdr3.responders

sum_total[Response == "Long-term Responder"]$Response <- "LtR"

sum_total[TimePoint_Final == "4 weeks after CAR T"]$TimePoint_Final <- "4wk post-CAR T"
sum_total[TimePoint_Final == "6 months after CAR T"]$TimePoint_Final <- "6mo post-CAR T"
sum_total[TimePoint_Final == "12 months after CAR T"]$TimePoint_Final <- "12mo post-CAR T"
sum_total$TimePoint_Final <- factor(sum_total$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))


p1 <- ggplot(data = sum_total, aes(x=TimePoint_Final, y=`Clonesize_proportion`, color=cdr3s_aa)) + 
  ggtitle("Cytotoxic Expanded TCRs 4wk post-CAR T (n=43)")+ ylab("Clonotypic Relative Abundance (%)")+ xlab("") + 
  geom_point(size=2, alpha=0.9, width = 0.05)+ 
   geom_line(aes(group=cdr3s_aa), linetype = "dashed", alpha = 0.5, colour = "black", data=sum_total) +
  scale_color_manual(values =  inferno(50)) + 
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,color="black",face="bold"), 
        axis.title=element_text(size=14,color="black",face="bold"),
        legend.position = 'none', 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.x = element_text(angle=65, hjust=1,size=12, face="bold"))
        
print(p1)

