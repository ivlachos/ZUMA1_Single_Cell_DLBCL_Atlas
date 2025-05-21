# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library("ggplot2")

# ===== Functions =====

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


# ===== Initialization =====

metadata <- readRDS("/data/cytotoxic_activity_scores_NKcells.txt")
setDT(metadata)

metadata$Response <- factor(metadata$Response, levels=c("Long-term Responder","Relapsed","non-Responder"))
metadata$TimePoint_Final <- factor(metadata$TimePoint_Final, levels=c("Leukapheresis","4 weeks after CAR T"))
metadata[TimePoint_Final == "Leukapheresis"]$TimePoint_Final <- "Leukapheresis"  
metadata[TimePoint_Final == "4 weeks after CAR T"]$TimePoint_Final <- "4wk post-CAR T"   


metadata$TimePoint_Final <- factor(metadata$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T"))

colors <- c("blue3","darkolivegreen3","brown3")

metadata <- metadata[TimePoint_Final %in% c("Leukapheresis","4wk post-CAR T")]
metadata$TimePoint_Final <- factor(metadata$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T"))


metadata[Response == "Long-term Responder"]$Response <- "LtR"
metadata[Response == "Relapsed"]$Response <- "R"
metadata[Response == "non-Responder"]$Response <- "NR"
metadata$Response <- factor(metadata$Response , levels = c("LtR","R","NR"))

# ===== Figure 4b =====

p <- ggplot(data = metadata[Response != "NR"], aes(x=`TimePoint_Final`, y=`cytotoxic_activity_UCell`, fill=Response)) + facet_wrap(~`Response`) +  
  ylab("Cytotoxic/Activation Signature")+ xlab("")+ 
  geom_boxplot(width=0.6, position=dodge,alpha=0.6,
               outlier.shape=NA)+scale_fill_manual(values = colors) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=12,color="black",face="bold"), 
        axis.title=element_text(size=12,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =12),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = 'top', 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=12, face="bold"),
        strip.text = element_text(size=12))


print(p)
