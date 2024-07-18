# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(networkdata)
library(igraph)
library(data.table)
library(ggraph)
library(graphlayouts)


# ===== Functions =====


plot_graph <- function(FC,limit_value,name){
  
  
  n = FC
  n <- n[,c("from", "to"), with=F]
  n <- n[to != "NA"]
  g <- graph_from_data_frame(n, directed=TRUE)
  
  ###########
  
  setnames(FC, "from", "name")
  FC$to <- NULL
  FC <- unique(FC)
  FC$Color <- "darkblue"
  FC[FC < 0]$Color <- "darkred"
  layout <- create_layout(g, layout="stress")
  layout <- merge(layout, FC, by="name")
  layout <- setDT(layout)
  layout <- layout[order(`.ggraph.index`)]
  
  g <- g %>% set_vertex_attr ("FC", value = layout$FC) %>% set_vertex_attr ("Color", value = layout$Color)
  
  p <- ggraph(g, layout = "kk") + ggtitle(name) + 
    geom_edge_link(edge_colour = "black", arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) +
    geom_node_point(aes(fill = FC),size=5,  shape = 21) +
    scale_fill_gradient2(low = "darkred", mid = "white",
                         high = "darkblue",name = "Log2(FC)",
                         limits=limit_value) + 
    
    geom_node_text(aes(label = name, color = Color), repel=TRUE, size=4) +
    
    scale_color_manual(values=c("darkblue","darkred")) + 
    
    scale_edge_width(range = c(0.2, 3)) +
    theme_graph() +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          legend.title = element_text(colour="black", size = 14, face = "bold"), title = element_text(colour="black", size = 18),
          legend.text = element_text(colour="black", size = 12),
          legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(0.5, 'cm'))
  
  p <- p + guides(color = "none")
  
  return(p)
}


# ===== Initialization =====

metadata <- readRDS("/data/metadata.RDS")


# ===== Trajectory plots =====

###### Remove samples with Unknown response and cells annotated as doublets. Retain only T Cells.

setDT(metadata)
metadata <- metadata[!metadata$orig.ident %in% c("R1","R2","T8"),]
metadata <- metadata[Response != "Unknown"]
metadata <- metadata[!Alias %in% c("DBL1","DBL2","DBL3","DBL4","DBL5","DBL6","DBL7","DBL8","DBL9")]

metadata$Label <- paste0(metadata$TimePoint_Final,"_",metadata$orig.ident,"@", metadata$Response)

metadata$Alias_initial <- metadata$Alias
countCells <- setnames(plyr::count(metadata$Alias), c("Alias", "No_of_Cells"))
metadata <- merge(metadata, countCells, by = "Alias")
annotation <- metadata
annotation <- unique(annotation[,c("Alias_initial","Alias","Compartment"), with=F])


metadata <- metadata[,c("Cell.id","orig.ident","TimePoint_Final","Patient_ID","cdr3s_aa","Label","Response", "Alias","Compartment"), with=F]
metadata <- metadata[Compartment %in% c("CD8 T Cells","CD4 T Cells","NK T Cells")]



sampleFreq <- setnames(setDT(plyr::count(metadata[,c("Label"), with=F])),
                       c("Label","All_Cells"))

metadata <-merge(metadata, sampleFreq, by = "Label")


########## Total Cellular proportions per timepoint 
total_proprotions_across_groups.1 <- setnames(setDT(plyr::count(metadata[,c("Label","TimePoint_Final","Alias"), with=F])),
                                              c("Label","TimePoint_Final","Alias","Counts"))

total_proprotions_across_groups <- dcast(total_proprotions_across_groups.1, Label+TimePoint_Final~Alias)

total_proprotions_across_groups$Condition <- gsub(".*_","", total_proprotions_across_groups$Label)
total_proprotions_across_groups[is.na(total_proprotions_across_groups)] <- 0

########## Calculate proportions
total_proprotions_across_groups.1 <- melt(total_proprotions_across_groups)
setnames(total_proprotions_across_groups.1, c("variable","value"), c("Alias","Counts"))

total_proprotions_across_groups.1 <- merge(total_proprotions_across_groups.1, sampleFreq, by="Label")

annotation <- annotation[,c("Alias","Compartment"), with=F]
total_proprotions_across_groups.1 <- merge(total_proprotions_across_groups.1, annotation, by="Alias")


total_proprotions_across_groups.1$Condition <- gsub(".*@","", total_proprotions_across_groups.1$Label)
total_proprotions_across_groups.1$Cellular_Proportion = (total_proprotions_across_groups.1$Counts/total_proprotions_across_groups.1$All_Cells)*100

setDT(total_proprotions_across_groups.1)
total_proprotions_across_groups.1 <- unique(total_proprotions_across_groups.1[,c("Alias","TimePoint_Final","Condition","Cellular_Proportion"), with=F])
total_proprotions_across_groups.1 <- total_proprotions_across_groups.1[,mean_CellProportions := mean(Cellular_Proportion), by = c("Alias,TimePoint_Final,Condition")]

total_proprotions_across_groups.1$Cellular_Proportion <- total_proprotions_across_groups.1$mean_CellProportions
total_proprotions_across_groups.1$mean_CellProportions <- NULL
total_proprotions_across_groups.1 <- unique(total_proprotions_across_groups.1)

total_proprotions_across_groups.1$Cellular_Proportion <- total_proprotions_across_groups.1$Cellular_Proportion + 1


total_proprotions_across_groups.1.leukapheresis <- total_proprotions_across_groups.1[TimePoint_Final == "Leukapheresis"]
total_proprotions_across_groups.1.4weeks <- total_proprotions_across_groups.1[TimePoint_Final == "4 weeks after CAR T"]
total_proprotions_across_groups.1.6months <- total_proprotions_across_groups.1[TimePoint_Final == "6 months after CAR T"]

total_proprotions_across_groups.1.leukapheresis <- dcast(total_proprotions_across_groups.1.leukapheresis, Alias~Condition)
total_proprotions_across_groups.1.4weeks <- dcast(total_proprotions_across_groups.1.4weeks, Alias~Condition)
total_proprotions_across_groups.1.6months <- dcast(total_proprotions_across_groups.1.6months, Alias~Condition)

total_proprotions_across_groups.1.leukapheresis$FC_Responders_vs_Relapsed <- total_proprotions_across_groups.1.leukapheresis$`Long-term Responder`/total_proprotions_across_groups.1.leukapheresis$`Relapsed`
total_proprotions_across_groups.1.leukapheresis$FC_Responders_vs_nonResponders <- total_proprotions_across_groups.1.leukapheresis$`Long-term Responder`/total_proprotions_across_groups.1.leukapheresis$`non-Responder`

total_proprotions_across_groups.1.4weeks$FC_Responders_vs_Relapsed <- total_proprotions_across_groups.1.4weeks$`Long-term Responder`/total_proprotions_across_groups.1.4weeks$`Relapsed`
total_proprotions_across_groups.1.4weeks$FC_Responders_vs_nonResponders <- total_proprotions_across_groups.1.4weeks$`Long-term Responder`/total_proprotions_across_groups.1.4weeks$`non-Responder`

total_proprotions_across_groups.1.6months$FC_Responders_vs_Relapsed <- total_proprotions_across_groups.1.6months$`Long-term Responder`/total_proprotions_across_groups.1.6months$`Relapsed`


total_proprotions_across_groups.1.6months$FC_Responders_vs_nonResponders <- 0

total_proprotions_across_groups.1.leukapheresis <- setDT(total_proprotions_across_groups.1.leukapheresis)[,c("Alias","FC_Responders_vs_Relapsed","FC_Responders_vs_nonResponders"), with=F]
total_proprotions_across_groups.1.leukapheresis$Condition <- "Leukapheresis"
total_proprotions_across_groups.1.4weeks <- setDT(total_proprotions_across_groups.1.4weeks)[,c("Alias","FC_Responders_vs_Relapsed","FC_Responders_vs_nonResponders"), with=F]
total_proprotions_across_groups.1.4weeks$Condition <- "4 weeks after CAR T"
total_proprotions_across_groups.1.6months <- setDT(total_proprotions_across_groups.1.6months)[,c("Alias","FC_Responders_vs_Relapsed","FC_Responders_vs_nonResponders"), with=F]
total_proprotions_across_groups.1.6months$Condition <- "6 months after CAR T"

total <- rbind(total_proprotions_across_groups.1.leukapheresis, total_proprotions_across_groups.1.4weeks,total_proprotions_across_groups.1.6months)


############################
########### CD8 Trajectories
############################


cd8 <- total[regexpr("CD8.*|NK T cells.*", total$Alias)==TRUE]

########### Fix the Graph
setnames(cd8, "Alias","from")

cd8$to <- c("NA","NA","CD8 Tem/MPECs","NK T cells/CD8+CD16+","NA","CD8 Tem/TGFb1+CXCR4+","CD8 Tem/Proinflammatory","NA","CD8 Tem/Highly cytotoxic","NA", "CD8 Ki67+/TOP2A+", "CD8 Tex/Ki67+BATF+","NK T cells/Proinflammatory","NA",
            "NA","NA","CD8 Tem/MPECs","NK T cells/CD8+CD16+","NA","CD8 Tem/TGFb1+CXCR4+","CD8 Tem/Proinflammatory","NA","CD8 Tem/Highly cytotoxic","NA", "CD8 Ki67+/TOP2A+", "CD8 Tex/Ki67+BATF+","NK T cells/Proinflammatory","NA",
            "NA","NA","CD8 Tem/MPECs","NK T cells/CD8+CD16+","NA","CD8 Tem/TGFb1+CXCR4+","CD8 Tem/Proinflammatory","NA","CD8 Tem/Highly cytotoxic","NA", "CD8 Ki67+/TOP2A+", "CD8 Tex/Ki67+BATF+","NK T cells/Proinflammatory","NA") 


cd8.new <- setnames(data.frame(c("CD8 Tem/TGFb1+CXCR4+","CD8 Tem/PF4+", "CD8 Tem/PF4+", "CD8 Tem/Highly cytotoxic", "CD8 Tem/Highly cytotoxic"),
                                 c("Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis"),
                                 c("CD8 Tem/PF4+","CD8 Temra/senescent","CD8 MAIT","CD8 Tem/LLEC","CD8 Tex/Ki67+IRF4+")),
                      c("from","Condition","to"))

cd8.new2 <- setnames(data.frame(c("CD8 Tem/TGFb1+CXCR4+","CD8 Tem/PF4+", "CD8 Tem/PF4+", "CD8 Tem/Highly cytotoxic", "CD8 Tem/Highly cytotoxic"),
                                  c("4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T"),
                                  c("CD8 Tem/PF4+","CD8 Temra/senescent","CD8 MAIT","CD8 Tem/LLEC","CD8 Tex/Ki67+IRF4+")),
                       c("from","Condition","to"))

cd8.new3 <- setnames(data.frame(c("CD8 Tem/TGFb1+CXCR4+","CD8 Tem/PF4+", "CD8 Tem/PF4+", "CD8 Tem/Highly cytotoxic", "CD8 Tem/Highly cytotoxic"),
                                  c("6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T"),
                                  c("CD8 Tem/PF4+","CD8 Temra/senescent","CD8 MAIT","CD8 Tem/LLEC","CD8 Tex/Ki67+IRF4+")),
                       c("from","Condition","to"))

cd8.initial <- cd8
cd8.initial$to <- NULL
cd8.new <- merge(cd8.new, cd8.initial, by = c("from","Condition"))
cd8.new2 <- merge(cd8.new2, cd8.initial, by = c("from","Condition"))
cd8.new3 <- merge(cd8.new3, cd8.initial, by = c("from","Condition"))

cd8.new <- rbind(cd8.new, cd8.new2,cd8.new3)
cd8.new <- setDT(cd8.new)
cd8.new <- cd8.new[,colnames(cd8), with=F]
cd8 <- rbind(cd8, cd8.new)

###### FC CD8 trajectories for Responder vs Relapsed

cd8$FC_Responders_vs_Relapsed <- log2(cd8$FC_Responders_vs_Relapsed)
cd8$FC_Responders_vs_nonResponders <- log2(cd8$FC_Responders_vs_nonResponders)

cd8.leukapheresis.relapsed <- cd8[Condition == "Leukapheresis"]
cd8.leukapheresis.relapsed$FC <- cd8.leukapheresis.relapsed$FC_Responders_vs_Relapsed
cd8.leukapheresis.relapsed <- plot_graph(cd8.leukapheresis.relapsed, c(min(cd8$FC_Responders_vs_Relapsed),max(cd8$FC_Responders_vs_Relapsed)), name="Leukapheresis")


cd8.4weeks.relapsed <- cd8[Condition == "4 weeks after CAR T"]
cd8.4weeks.relapsed$FC <- cd8.4weeks.relapsed$FC_Responders_vs_Relapsed
cd8.4weeks.relapsed <- plot_graph(cd8.4weeks.relapsed,c(min(cd8$FC_Responders_vs_Relapsed), max(cd8$FC_Responders_vs_Relapsed)), "4 weeks after CAR T")


cd8.6months.relapsed <- cd8[Condition == "6 months after CAR T"]
cd8.6months.relapsed$FC <- cd8.6months.relapsed$FC_Responders_vs_Relapsed
cd8.6months.relapsed <- plot_graph(cd8.6months.relapsed,c(min(cd8$FC_Responders_vs_Relapsed),max(cd8$FC_Responders_vs_Relapsed)),"6 months after CAR T")


####### Plot the CD8 trajectories per timepoint for Responder vs Relapsed


print(cd8.leukapheresis.relapsed)

print(cd8.4weeks.relapsed)

print(cd8.6months.relapsed)


###### FC CD8 trajectories for Responder vs non-Responder


cd8.leukapheresis.Nonresponder <- cd8[Condition == "Leukapheresis"]
cd8.leukapheresis.Nonresponder$FC <- cd8.leukapheresis.Nonresponder$FC_Responders_vs_nonResponders
cd8.leukapheresis.Nonresponder <- plot_graph(cd8.leukapheresis.Nonresponder,c(min(cd8[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders),
                                                                              max(cd8[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders)), "Leukapheresis")


cd8.4weeks.Nonresponder <- cd8[Condition == "4 weeks after CAR T"]
cd8.4weeks.Nonresponder$FC <- cd8.4weeks.Nonresponder$FC_Responders_vs_nonResponders
cd8.4weeks.Nonresponder <- plot_graph(cd8.4weeks.Nonresponder,c(min(cd8[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders),
                                                                max(cd8[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders)), "4 weeks after CAR T")


####### Plot the CD8 trajectories per timepoint for Responder vs non-Responder
print(cd8.leukapheresis.Nonresponder)

print(cd8.4weeks.Nonresponder)




############################
########### CD4 Trajectories
############################


cd4 <- total[regexpr("CD4.*|T regs", total$Alias)==TRUE]


########### Fix the Graph
setnames(cd4, "Alias","from")

cd4$to <- c("NA","CD4 Teh/TCF7+LEF1+","NA","CD4 Th2","NA","NA","NA","T regs","CD4 Th17/T regs","NA","NA","CD4 Th1","CD4 Th17","NA",
            "NA","CD4 Teh/TCF7+LEF1+","NA","CD4 Th2","NA","NA","NA","T regs","CD4 Th17/T regs","NA","NA","CD4 Th1","CD4 Th17","NA",
            "NA","CD4 Teh/TCF7+LEF1+","NA","CD4 Th2","NA","NA","NA","T regs","CD4 Th17/T regs","NA","NA","CD4 Th1","CD4 Th17","NA") 


cd4.new <- setnames(data.frame(c("CD4 Th2/Th17","CD4 Teh/TCF7+LEF1+","CD4 Th1","CD4 Th1","CD4 Th1","CD4 Th2","CD4 Th2"),
                                 c("Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis","Leukapheresis"),
                                 c("CD4 Th17/BCL6+","CD4 Th2/Th17","CD4 Teh/TNF+","CD4 CTL","CD4 Tex/Ki67+","CD4 Tem/Teh","CD4 Teh/PF4+")),
                      c("from","Condition","to"))

cd4.new2 <- setnames(data.frame(c("CD4 Th2/Th17","CD4 Teh/TCF7+LEF1+","CD4 Th1","CD4 Th1","CD4 Th1","CD4 Th2","CD4 Th2"),
                                  c("4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T","4 weeks after CAR T"),
                                  c("CD4 Th17/BCL6+","CD4 Th2/Th17","CD4 Teh/TNF+","CD4 CTL","CD4 Tex/Ki67+","CD4 Tem/Teh","CD4 Teh/PF4+")),
                       c("from","Condition","to"))

cd4.new3 <- setnames(data.frame(c("CD4 Th2/Th17","CD4 Teh/TCF7+LEF1+","CD4 Th1","CD4 Th1","CD4 Th1","CD4 Th2","CD4 Th2"),
                                  c("6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T","6 months after CAR T"),
                                  c("CD4 Th17/BCL6+","CD4 Th2/Th17","CD4 Teh/TNF+","CD4 CTL","CD4 Tex/Ki67+","CD4 Tem/Teh","CD4 Teh/PF4+")),
                       c("from","Condition","to"))


cd4.initial <- cd4
cd4.initial$to <- NULL
cd4.new <- merge(cd4.new, cd4.initial, by = c("from","Condition"))
cd4.new2 <- merge(cd4.new2, cd4.initial, by = c("from","Condition"))
cd4.new3 <- merge(cd4.new3, cd4.initial, by = c("from","Condition"))

cd4.new <- rbind(cd4.new, cd4.new2,cd4.new3)
cd4.new <- setDT(cd4.new)
cd4.new <- cd4.new[,colnames(cd4), with=F]
cd4 <- rbind(cd4, cd4.new)

cd4$FC_Responders_vs_Relapsed <- log2(cd4$FC_Responders_vs_Relapsed)
cd4$FC_Responders_vs_nonResponders <- log2(cd4$FC_Responders_vs_nonResponders)

cd4.leukapheresis.relapsed <- cd4[Condition == "Leukapheresis"]
cd4.leukapheresis.relapsed$FC <- cd4.leukapheresis.relapsed$FC_Responders_vs_Relapsed
cd4.leukapheresis.relapsed <- plot_graph(cd4.leukapheresis.relapsed,limit_value =c(min(cd4$FC_Responders_vs_Relapsed),max(cd4$FC_Responders_vs_Relapsed)), name="Leukapheresis")


cd4.4weeks.relapsed <- cd4[Condition == "4 weeks after CAR T"]
cd4.4weeks.relapsed$FC <- cd4.4weeks.relapsed$FC_Responders_vs_Relapsed
cd4.4weeks.relapsed <- plot_graph(cd4.4weeks.relapsed,c(min(cd4$FC_Responders_vs_Relapsed),max(cd4$FC_Responders_vs_Relapsed)), "4 weeks after CAR T")


cd4.6months.relapsed <- cd4[Condition == "6 months after CAR T"]
cd4.6months.relapsed$FC <- cd4.6months.relapsed$FC_Responders_vs_Relapsed
cd4.6months.relapsed <- plot_graph(cd4.6months.relapsed, c(min(cd4$FC_Responders_vs_Relapsed),max(cd4$FC_Responders_vs_Relapsed)),"6 months after CAR T")

####### Plot the CD4 trajectories per timepoint for Responder vs Relapsed


print(cd4.leukapheresis.relapsed)
print(cd4.4weeks.relapsed)
print(cd4.6months.relapsed)


####### FC Responder vs non Responder

cd4.leukapheresis.Nonresponder <- cd4[Condition == "Leukapheresis"]
cd4.leukapheresis.Nonresponder$FC <- cd4.leukapheresis.Nonresponder$FC_Responders_vs_nonResponders
cd4.leukapheresis.Nonresponder <- plot_graph(cd4.leukapheresis.Nonresponder,c(min(cd4[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders),
                                                                                       max(cd4[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders)), "Leukapheresis")


cd4.4weeks.Nonresponder <- cd4[Condition == "4 weeks after CAR T"]
cd4.4weeks.Nonresponder$FC <- cd4.4weeks.Nonresponder$FC_Responders_vs_nonResponders
cd4.4weeks.Nonresponder <- plot_graph(cd4.4weeks.Nonresponder,c(min(cd4[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders),
                                                                         max(cd4[Condition == "Leukapheresis"|Condition == "4 weeks after CAR T"]$FC_Responders_vs_nonResponders)), "4 weeks after CAR T")


####### Plot the cd4 trajectories per timepoint for Responder vs non-Responder


print(cd4.leukapheresis.Nonresponder)

print(cd4.4weeks.Nonresponder)
