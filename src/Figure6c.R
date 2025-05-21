# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(pROC)
library(ggplot2)
library(grid)


# ===== Figure 6c =====
roclist <- readRDS("/data/Figure6c.RDS")

rocplot <- ggroc(roclist, aes='colour', legacy.axes = F, aes(size=Folds)) +
    scale_size_manual(values = c(1,1,1,1)) +
    geom_segment(aes(x=1, xend=0, y=0, yend=1), color='grey2', linetype='dashed') +
  scale_color_manual(values = c("blue","orange","purple","brown")) +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(size=10,color="black",face="bold"), 
          axis.title=element_text(size=10,color="black",face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =11),
          legend.key = element_rect(fill = "transparent", linewidth = 1,colour = "transparent"),
          legend.key.size = unit(0.5, "cm"), 
          # legend.position = 'right', 
          legend.position = c(0.6,0.12), 
          axis.line = element_line(size = 0.5, colour = "black"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=10, face="bold"),
          axis.text.x = element_text(size=11, face="bold"))
    
  
print(rocplot)
  