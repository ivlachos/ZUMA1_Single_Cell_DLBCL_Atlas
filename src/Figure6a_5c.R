# R version 3.6.0 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(ggplot2)
library(dplyr)
library(ggbreak)
library(plyr)

# ===== Initialization =====


clonotypes_with_Type_4wk <- fread("/data/TCRs_Leukapheresis_4 weeks_CytotoxicCombined.txt")
clonotypes_with_Type_4wk <- clonotypes_with_Type_4wk[,c("Leukapheresis","4 weeks after CAR T"):=NULL]
clonotypes_with_Type_4wk$Timepoint <- "4wk post-CAR T"


clonotypes_with_Type_6mo <- fread("/data/TCRs_Leukapheresis_6 months_CytotoxicCombined.txt")
clonotypes_with_Type_6mo <- clonotypes_with_Type_6mo[,c("Leukapheresis","6 months after CAR T"):=NULL]
clonotypes_with_Type_6mo$Timepoint <- "6mo post-CAR T"


clonotypes_with_Type_12mo <- fread("/data/TCRs_Leukapheresis_12 months_CytotoxicCombined.txt")
clonotypes_with_Type_12mo <- clonotypes_with_Type_12mo[,c("Leukapheresis","12 months after CAR T"):=NULL]
clonotypes_with_Type_12mo$Timepoint <- "12mo post-CAR T"

clonotypes_with_Type <- rbind(clonotypes_with_Type_4wk, clonotypes_with_Type_6mo, clonotypes_with_Type_12mo)

clonotypes_with_Type <- clonotypes_with_Type[Response == "Long-term Responder"]
clonotypes_with_Type[Response == "Long-term Responder"]$Response <- "LtR"


# ===== Figure 6a =====

clonotypes_with_Type_selected <- clonotypes_with_Type[Alias %in% c("T cells cytotoxic")]
data <- clonotypes_with_Type_selected

shape_values <- vector()
plot_labels <- vector()
color_values <- vector()

data[Type %in% c("Persistent","Contracted")]$Type <- "Persistent/Contracted"

if(sum(data$Type == "Persistent/Contracted") > 0){
  shape_values <- c(3)
  plot_labels <- c("Persistent/Contracted")
  color_values <- c("gray48")
  
}

if(sum(data$Type == "Novel") > 0){
  shape_values <- c(shape_values, 15)
  plot_labels <- c(plot_labels, "Novel")
  color_values <- c(color_values,"#21908CFF")
  
}
if(sum(data$Type == "Expanded") > 0){
  shape_values <- c(shape_values, 17)
  plot_labels <- c(plot_labels, "Expanded")
  color_values <- c(color_values,"brown")
  
}

data$Type <- factor(data$Type, levels = c("Persistent/Contracted","Novel","Expanded"))


data$Timepoint <- factor(data$Timepoint, levels = c("Leukapheresis","4wk post-CAR T","6mo post-CAR T","12mo post-CAR T"))
data$Response <- factor(data$Response, levels = c("LtR","R","NR"))


library(ggbreak)

data$Pre_frequency <- data$Pre_frequency*100
data$Post_frequency <- data$Post_frequency*100

p <- ggplot(data, aes(x=Pre_frequency, y=Post_frequency, color = Type, shape = Type)) +   facet_grid(rows = vars(Response), cols = vars(Timepoint), scales = "free_x") + 
  geom_point(data = subset(data, Type != "Expanded"), size = 3, alpha = 1) +
  geom_point(data = subset(data, Type == "Expanded"), size = 3, alpha =1) + 
  scale_color_manual(values = color_values, labels = plot_labels) + 
  scale_shape_manual(values = c(3,15, 17)) + 
  ylab(paste0("post-CAR T, Clonosize proportion (%)")) +

  xlab(paste0("Leukapheresis, Clonosize proportion (%)")) + theme( panel.background = element_rect(fill = "white"),
                                                                   panel.grid.major = element_line(color = "grey80"),
                                                                   panel.grid.minor = element_line(color = "grey90"),     
                                                                  plot.title = element_text(size=15,color="black",face="bold"), 
                                                                  axis.title=element_text(size=14,color="black",face="bold"),
                                                                  legend.title = element_blank(),
                                                                  axis.text.x = element_text(size=12),axis.text.y = element_text(size=14),
                                                                  legend.text = element_text(colour="black", size = 15),
                                                                  legend.key = element_rect(fill = "transparent", size = 15,colour = "transparent"),
                                                                  legend.key.size = unit(0.2, "cm"), #legend.position = "none",
                                                                  legend.position = "top", 
                                                                  legend.key.width = unit(1,"cm"), 
                                                                  axis.line = element_line(size = 0.5, colour = "black"),
                                                                  plot.margin = unit(c(0.1,0.7,0.1,0.1), "cm"), strip.text = element_text(size=13)) 

print(p)


# ===== Figure 5c =====


clonotypes_with_Type <- rbind(clonotypes_with_Type_4wk, clonotypes_with_Type_6mo, clonotypes_with_Type_12mo)

clonotypes_with_Type[Response == "Long-term Responder"]$Response <- "LtR"
clonotypes_with_Type[Response == "non-Responder"]$Response <- "NR"
clonotypes_with_Type[Response == "Relapsed"]$Response <- "R"

clonotypes_with_Type_selected <- clonotypes_with_Type[Alias %in% c("T cells cytotoxic")]


data <- clonotypes_with_Type_selected

data[Type %in% c("Persistent","Contracted")]$Type <- "Persistent/Contracted"

data <- data[!((Timepoint == "12mo post-CAR T" & Response == "R") | (Timepoint == "6mo post-CAR T" & Response == "R"))]

data$Timepoint <- factor(data$Timepoint, levels = c("Leukapheresis","4wk post-CAR T","6mo post-CAR T","12mo post-CAR T"))
data$Response <- factor(data$Response, levels = c("LtR","R","NR"))

sum_cdr3 <- plyr::count(data[,c("Response","Type","Timepoint","Patient")])
setDT(sum_cdr3)
sum_cdr3.dc <- dcast(sum_cdr3,Response+Timepoint+Patient~Type)
setDT(sum_cdr3.dc)
sum_cdr3.dc[is.na(sum_cdr3.dc$`Expanded`)]$`Expanded` <- 0
sum_cdr3.dc <- melt(sum_cdr3.dc)
setDT(sum_cdr3.dc)
sum_cdr3.dc <- sum_cdr3.dc[variable == "Expanded"]

sum_cdr3.dc$freq <- 0
sum_cdr3.dc[value > 1]$freq <- 1
sum_cdr3.dc <- sum_cdr3.dc[Timepoint == "4wk post-CAR T"]
# Filter the rows where freq = 1
sum_cdr3.dc <- sum_cdr3.dc[, .(count = sum(freq == 1), total = .N), by = Response]
# Calculate the percentage
sum_cdr3.dc[, percentage := (count / total) * 100]

colors <- c("blue3","darkolivegreen3","#d95f0e")

sum_cdr3.dc$Response <- factor(sum_cdr3.dc$Response, levels = c("LtR","R","NR"))

sum_cdr3.dc$Timepoint <- "4wk post-CAR T"

plot <- ggplot(sum_cdr3.dc, aes(y = percentage, x = Response, fill = Response)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge(0.9), alpha = 0.6, width = 0.5) +
  ylab("Patients with Expanded Cytotoxic TCRs (%)") +scale_fill_manual(values = colors) + xlab("") + 
  theme_minimal() + ggtitle("4wk post-CAR T") +
  theme( panel.background = element_rect(fill = "white"),
         panel.grid.minor = element_line(color = "grey90"),     
         plot.title = element_text(size=9,color="black",face="bold"), 
         axis.title=element_text(size=9,color="black",face="bold"),
         legend.title = element_blank(),
         axis.text.x = element_blank(),
         legend.text = element_text(colour="black", size = 9),
         legend.key = element_rect(fill = "transparent", size = 9,colour = "transparent"),
         legend.key.size = unit(0.1, "cm"), #legend.position = "none",
         legend.position = "right", 
         legend.key.width = unit(0.5,"cm"), 
         axis.line = element_line(size = 0.5, colour = "black"),
         plot.margin = unit(c(0.1,0.7,0.1,0.1), "cm"), strip.text = element_text(size=11)) 


print(plot)

