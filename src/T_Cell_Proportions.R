# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(dplyr)
library(viridis)
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(reshape2)
library(plyr)

# ===== Functions =====

data_summary <- function(data, varname="Clonosize", 
                         groupnames=c("Response","Alias")){
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


# ===== Initialization =====

metadata <- readRDS("/data/metadata.RDS")


###### Remove samples with Unknown response and cells annotated as doublets. Retain only T Cells.

setDT(metadata)
metadata <- metadata[!metadata$orig.ident %in% c("R1","R2","T8"),]
metadata <- metadata[Response != "Unknown"]
metadata <- metadata[!Alias %in% c("DBL1","DBL2","DBL3","DBL4","DBL5","DBL6","DBL7","DBL8","DBL9")]
metadata <- metadata[!(Patient_ID == "P10" & TimePoint_Final == "12mo post-CAR T")]
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


###### Calculate proportions

total_proprotions_across_groups.1 <- setnames(setDT(plyr::count(metadata[,c("Label","TimePoint_Final","Alias"), with=F])),
                                              c("Label","TimePoint_Final","Alias","Counts"))

total_proprotions_across_groups <- dcast(total_proprotions_across_groups.1, Label+TimePoint_Final~Alias)

total_proprotions_across_groups$Response <- gsub(".*@","", total_proprotions_across_groups$Label)
total_proprotions_across_groups[is.na(total_proprotions_across_groups)] <- 0

total_proprotions_across_groups.1 <- melt(total_proprotions_across_groups)
setnames(total_proprotions_across_groups.1, c("variable","value"), c("Alias","Counts"))

total_proprotions_across_groups.1 <- merge(total_proprotions_across_groups.1, sampleFreq, by="Label")

annotation <- annotation[,c("Alias","Compartment"), with=F]
total_proprotions_across_groups.1 <- merge(total_proprotions_across_groups.1, annotation, by="Alias")

total_proprotions_across_groups.1$Response <- gsub(".*@","", total_proprotions_across_groups.1$Label)
total_proprotions_across_groups.1$Cellular_Proportion = (total_proprotions_across_groups.1$Counts/total_proprotions_across_groups.1$All_Cells)*100

total_proprotions_across_groups <- setDT(total_proprotions_across_groups)


# ===== CD4 cell proportions Responders vs R=====
setDT(total_proprotions_across_groups.1)
data <- total_proprotions_across_groups.1[regexpr("CD4.*|T regs",total_proprotions_across_groups.1$Alias)==TRUE]
data$Alias <- factor(data$Alias)
data <- data[Response != "NR"]


data$Response <- factor(data$Response, levels = c("LtR", "R"))
df2 <- data_summary(data, varname=c("Cellular_Proportion"), 
                    groupnames=c("Response","TimePoint_Final","Alias"))


data$TimePoint_Final <- factor(data$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "R"))

df2$Alias <- factor(df2$Alias, levels = c("CD4 CTL","T regs","CD4 Tcm/naive","CD4 Teh/TCF7+LEF1+","CD4 Tem/Teh","CD4 Th1","CD4 Th2","CD4 Th2/Th17",
                                          "CD4 Th17","CD4 Th17/BCL6+","CD4 Th17/T regs","CD4 Tex/Ki67+","CD4 Teh/PF4+","CD4 Teh/TNF+"))

setDT(df2)

data.plot <-  ggplot(data = df2, aes(x=TimePoint_Final, y=`Cellular_Proportion`, fill=Response)) + facet_wrap(~Alias, ncol=7, scales = "free_y")+
  geom_bar(stat="identity", color="black",position=position_dodge(0.9), alpha=0.6, width=0.8) + scale_alpha_manual(values = c("0.5"=0.4, "1"=1), guide='none')+
  geom_errorbar(aes(ymin=Cellular_Proportion-se, ymax=Cellular_Proportion+se), width=.2,
                position=position_dodge(.9)) +  
  geom_jitter(data = data, aes( x = TimePoint_Final, y =Cellular_Proportion,
                                fill =factor(Response, levels = c( "LtR","R"))),
              position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=1) + 
  ggtitle("")+ ylab("Cellular Proportion (%)")+ xlab("") + 
  scale_fill_manual(values = c("blue3","darkolivegreen3")) + scale_y_continuous(expand = c(0, 0), limits = function(x) c(min(x), max(x) + 3))+
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17,color="black",face="bold"), 
        axis.title=element_text(size=18,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =17),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = "top", 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=17, face="bold"),
        strip.text = element_text(size=14, face="bold"),strip.background = element_rect(colour="black", fill="white"))


print(data.plot)


# ===== CD4 cell proportions Responders vs NR=====
setDT(total_proprotions_across_groups.1)
data <- total_proprotions_across_groups.1[regexpr("CD4.*|T regs",total_proprotions_across_groups.1$Alias)==TRUE]
data$Alias <- factor(data$Alias)
data <- data[Response != "R"]
data <- data[TimePoint_Final %in% c("Leukapheresis","4wk post-CAR T")]


data$Response <- factor(data$Response, levels = c("LtR", "NR"))
df2 <- data_summary(data, varname=c("Cellular_Proportion"), 
                    groupnames=c("Response","TimePoint_Final","Alias"))


data$TimePoint_Final <- factor(data$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "NR"))

df2$Alias <- factor(df2$Alias, levels = c("CD4 CTL","T regs","CD4 Tcm/naive","CD4 Teh/TCF7+LEF1+","CD4 Tem/Teh","CD4 Th1","CD4 Th2","CD4 Th2/Th17",
                                          "CD4 Th17","CD4 Th17/BCL6+","CD4 Th17/T regs","CD4 Tex/Ki67+","CD4 Teh/PF4+","CD4 Teh/TNF+"))

setDT(df2)

data.plot <-  ggplot(data = df2, aes(x=TimePoint_Final, y=`Cellular_Proportion`, fill=Response)) + facet_wrap(~Alias, ncol=7, scales = "free_y")+
  geom_bar(stat="identity", color="black",position=position_dodge(0.9), alpha=0.6, width=0.8) + scale_alpha_manual(values = c("0.5"=0.4, "1"=1), guide='none')+
  geom_errorbar(aes(ymin=Cellular_Proportion-se, ymax=Cellular_Proportion+se), width=.2,
                position=position_dodge(.9)) +  
  geom_jitter(data = data, aes( x = TimePoint_Final, y =Cellular_Proportion,
                                fill =factor(Response, levels = c( "LtR","NR"))),
              position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=1) + 
  ggtitle("")+ ylab("Cellular Proportion (%)")+ xlab("") + 
  scale_fill_manual(values =  c("blue3","brown")) + 
  scale_y_continuous(expand = c(0, 0), limits = function(x) c(min(x), max(x) + 3))+
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17,color="black",face="bold"), 
        axis.title=element_text(size=18,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =17),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = "top", 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=17, face="bold"),
        strip.text = element_text(size=14, face="bold"),strip.background = element_rect(colour="black", fill="white"))


print(data.plot)


# ===== CD8 cell proportions Responders vs NR=====
setDT(total_proprotions_across_groups.1)
data <- total_proprotions_across_groups.1[regexpr("CD8.*|NK T.*",total_proprotions_across_groups.1$Alias)==TRUE]
data$Alias <- factor(data$Alias)
data <- data[Response != "NR"]


data$Response <- factor(data$Response, levels = c("LtR", "R"))
df2 <- data_summary(data, varname=c("Cellular_Proportion"), 
                    groupnames=c("Response","TimePoint_Final","Alias"))


data$TimePoint_Final <- factor(data$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "R"))


df2$Alias <- factor(df2$Alias, levels = c("CD8 Tem/Highly cytotoxic","CD8 Tem/TGFb1+CXCR4+", "CD8 Tem/Proinflammatory", "CD8 Ki67+/TOP2A+","CD8 Tex/Ki67+IRF4+","CD8 Tex/Ki67+BATF+",
                                          "CD8 Tcm/naive", "CD8 Tem/MPECs", "CD8 Temra/senescent","CD8 MAIT","CD8 Tem/LLEC","CD8 Tem/PF4+",
                                          "NK T cells/CD8+CD16+","NK T cells/Proinflammatory"))

setDF(df2)

data.plot <- ggplot(data = df2, aes(x=TimePoint_Final, y=`Cellular_Proportion`, fill=Response)) + facet_wrap(~Alias, ncol=4, scales = "free_y")+
  geom_bar(stat="identity", color="black",position=position_dodge(0.9), alpha=0.6, width=0.8) + scale_alpha_manual(values = c("0.5"=0.4, "1"=1), guide='none')+
  geom_errorbar(aes(ymin=Cellular_Proportion-se, ymax=Cellular_Proportion+se), width=.2,
                position=position_dodge(.9)) +  
  geom_jitter(data = data, aes( x = TimePoint_Final, y =Cellular_Proportion,
                                fill =factor(Response, levels = c( "LtR","R"))),
              position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=1) + 
  ggtitle("")+ ylab("Cellular Proportion (%)")+ xlab("") + 
  scale_fill_manual(values = c("blue3","darkolivegreen3")) + scale_y_continuous(expand = c(0, 0), limits = function(x) c(min(x), max(x) + 3))+
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17,color="black",face="bold"), 
        axis.title=element_text(size=18,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =17),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = "top", 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=17, face="bold"),
        strip.text = element_text(size=14, face="bold"),strip.background = element_rect(colour="black", fill="white"))


print(data.plot)



# ===== CD8 cell proportions Responders vs R=====
setDT(total_proprotions_across_groups.1)
data <- total_proprotions_across_groups.1[regexpr("CD8.*|NK T.*",total_proprotions_across_groups.1$Alias)==TRUE]
data$Alias <- factor(data$Alias)
data <- data[Response != "R"]
data <- data[TimePoint_Final %in% c("Leukapheresis","4wk post-CAR T")]


data$Response <- factor(data$Response, levels = c("LtR", "NR"))
df2 <- data_summary(data, varname=c("Cellular_Proportion"), 
                    groupnames=c("Response","TimePoint_Final","Alias"))


data$TimePoint_Final <- factor(data$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T", "6mo post-CAR T", "12mo post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "NR"))


df2$Alias <- factor(df2$Alias, levels = c("CD8 Tem/Highly cytotoxic","CD8 Tem/TGFb1+CXCR4+", "CD8 Tem/Proinflammatory", "CD8 Ki67+/TOP2A+","CD8 Tex/Ki67+IRF4+","CD8 Tex/Ki67+BATF+",
                                          "CD8 Tcm/naive", "CD8 Tem/MPECs", "CD8 Temra/senescent","CD8 MAIT","CD8 Tem/LLEC","CD8 Tem/PF4+",
                                          "NK T cells/CD8+CD16+","NK T cells/Proinflammatory"))

setDF(df2)

data.plot <- ggplot(data = df2, aes(x=TimePoint_Final, y=`Cellular_Proportion`, fill=Response)) + facet_wrap(~Alias, ncol=7, scales = "free_y")+
  geom_bar(stat="identity", color="black",position=position_dodge(0.9), alpha=0.6, width=0.8) + scale_alpha_manual(values = c("0.5"=0.4, "1"=1), guide='none')+
  geom_errorbar(aes(ymin=Cellular_Proportion-se, ymax=Cellular_Proportion+se), width=.2,
                position=position_dodge(.9)) +  
  geom_jitter(data = data, aes( x = TimePoint_Final, y =Cellular_Proportion,
                                fill =factor(Response, levels = c( "LtR","R"))),
              position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=1) + 
  ggtitle("")+ ylab("Cellular Proportion (%)")+ xlab("") + 
  scale_fill_manual(values = c("blue3","brown")) + scale_y_continuous(expand = c(0, 0), limits = function(x) c(min(x), max(x) + 3))+
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17,color="black",face="bold"), 
        axis.title=element_text(size=18,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =17),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = "top", 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=17, face="bold"),
        strip.text = element_text(size=14, face="bold"),strip.background = element_rect(colour="black", fill="white"))


print(data.plot)