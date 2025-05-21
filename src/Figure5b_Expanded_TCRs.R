# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(ggplot2)
library(plyr)

# ===== Functions =====

data_summary <- function(data, varname, groupnames){
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


# ===== Figure 5b =====

###### Remove samples with Unknown response and cells annotated as doublets. Retain only T Cells.

setDT(metadata)
metadata <- metadata[!metadata$orig.ident %in% c("R1","R2","T8"),]
metadata <- metadata[Response != "Unknown"]
metadata <- metadata[!Alias %in% c("DBL1","DBL2","DBL3","DBL4","DBL5","DBL6","DBL7","DBL8","DBL9")]
metadata <- metadata[!(Patient_ID == "P10" & TimePoint_Final == "12mo post-CAR T") & TimePoint_Final != "6mo post-CAR T" & TimePoint_Final != "12mo post-CAR T"]
metadata$Label <- paste0(metadata$TimePoint_Final,"_",metadata$orig.ident,"@", metadata$Response)
metadata <- metadata[Compartment %in% c("CD8 T Cells","CD4 T Cells","NK T Cells") & !is.na(cdr3s_aa)]

metadata <- metadata[,c("Cell.id","orig.ident","TimePoint_Final","Patient_ID","cdr3s_aa","Label","Response", "Alias","Compartment"), with=F]


sampleFreq <- setnames(setDT(plyr::count(metadata[,c("Label"), with=F])),
                       c("Label","All_Cells"))

metadata <-merge(metadata, sampleFreq, by = "Label")

annotation <- metadata
annotation <- unique(annotation[,c("Alias","Compartment"), with=F])
metadata$Response <- factor(metadata$Response, c("NR","R","LtR"))


######## Number of Expanded Clonotypes with TCR freq > 1

expanded_TCRs <- setnames(setDT(plyr::count(metadata[,c("Label","TimePoint_Final","Alias","cdr3s_aa"), with=F])),
                                                         c("Label","TimePoint_Final","Alias","cdr3s_aa","Counts"))

expanded_TCRs <- expanded_TCRs[Counts > 1]

expanded_TCRs <- expanded_TCRs[, Counts := sum(Counts), by=c("Label","Alias")]
expanded_TCRs$cdr3s_aa <- NULL
expanded_TCRs <- unique(expanded_TCRs)

expanded_TCRs <- dcast(expanded_TCRs, Label+TimePoint_Final~Alias)

expanded_TCRs$Response <- gsub(".*_","", expanded_TCRs$Label)
expanded_TCRs[is.na(expanded_TCRs)] <- 0

expanded_TCRs <- melt(expanded_TCRs)
setnames(expanded_TCRs, c("variable","value"), c("Alias","Counts"))

expanded_TCRs <- merge(expanded_TCRs, sampleFreq, by="Label")

annotation <- annotation[,c("Alias","Compartment"), with=F]
expanded_TCRs <- merge(expanded_TCRs, annotation, by="Alias")


expanded_TCRs$Response <- gsub(".*@","", expanded_TCRs$Label)
expanded_TCRs$ExpandedTCRs = (expanded_TCRs$Counts/expanded_TCRs$All_Cells)*100

###############################

expanded_TCRs <- setDT(expanded_TCRs)
data <- expanded_TCRs[Alias %in% c("CD8 Tem/Highly cytotoxic", "CD8 Tem/Proinflammatory", "CD8 Tem/TGFb1+CXCR4+", "CD8 Ki67+/TOP2A+","CD8 Tex/Ki67+IRF4+","CD4 CTL")]
data$Alias <- factor(data$Alias)

data$Response <- factor(data$Response, levels = c("LtR", "R","NR"))

df2 <- data_summary(data, varname="ExpandedTCRs",
                    groupnames=c("Response","TimePoint_Final","Alias"))

df2$Response <- factor(df2$Response, levels = c("LtR", "R","NR"))

if(length(unique(df2$Response))>2){
  colors <- c("blue3","darkolivegreen3","brown3")
}else{
  colors <- c("blue3","darkolivegreen3")
}

df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "R","NR"))
setDT(df2)


df2$TimePoint_Final <- factor(df2$TimePoint_Final, levels = c("Leukapheresis","4wk post-CAR T"))
df2$Response <- factor(df2$Response, levels = c("LtR", "R","NR"))

Tcell.plot <- ggplot(data = df2, aes(x=TimePoint_Final, y=`ExpandedTCRs`, fill=Response)) + facet_wrap(~Alias, ncol=2, scales = "free_y")+
  geom_bar(stat="identity", color="black",position=position_dodge(0.9), alpha=0.6, width=0.5) + scale_alpha_manual(values = c("0.5"=0.4, "1"=1), guide='none')+
  geom_errorbar(aes(ymin=ExpandedTCRs-se, ymax=ExpandedTCRs+se), width=.2,
                position=position_dodge(.9)) +  
  geom_jitter(data = data, aes( x = TimePoint_Final, y =ExpandedTCRs,
                                  fill =factor(Response, levels = c( "LtR","R","NR"))),
              position = position_dodge(width = 0.9), stat = "identity",alpha=0.4, size=1) + 
  ggtitle("")+ ylab("Expanded TCRs (%)")+ xlab("") + 
  scale_fill_manual(values = colors) + 
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17,color="black",face="bold"), 
        axis.title=element_text(size=18,color="black",face="bold"),
        legend.title =element_blank(),
        legend.text = element_text(colour="black", size =16),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = "top", 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 65, hjust = 1,size=17, face="bold"),
        strip.text = element_text(size=13, face="bold"),strip.background = element_rect(colour="black", fill="white"))+ guides(fill=guide_legend(title="TCR-clonotype freq > 4"))


print(Tcell.plot )

