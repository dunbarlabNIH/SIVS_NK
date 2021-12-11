## Script to make all figures for SIVS NK  manuscript
# Ryland Mortlock
# October 2021

## Install required libraries ###################################################
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("dplyr")
# install.packages("ggpubr")
# install.packages("RColorBrewer")
# install.packages("networkD3")
# install.packages("plyr")
# install.packages("pheatmap")

## Load required libraries ###################################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(networkD3)
library(plyr)
library(pheatmap)

# Set directories
data_dir <- "data_files_for_manuscript"
results_dir <- "R_plots"

# Function to get summary statistic
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary_sd <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

data_summary_se <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(!is.na(x[[col]]))))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}




## Box plot of NK proportions across tissues ###################################################
# Load subset proportionos data
NK_prop <- read.csv(file =  file.path(data_dir, "NK_proportion_of_lymphocytes.csv"))

# Set AxLN tissues as LN for SIVS SAMPLES
NK_prop$Tissue[grepl("Ax_LN",NK_prop$Tissue)] <- "Lymph Node"

# Order factors
NK_prop$Tissue <- factor(NK_prop$Tissue, levels = c("Blood","Bone Marrow","Lymph Node","Spleen",
                                                    "Liver","Lung", "Jejunum",
                                                    "BAL", "Colon", "Ing_LN", "Mes_LN"))

# Only keep the tissues we want to show
NK_prop <- NK_prop[NK_prop$Tissue %in% c("Blood","Bone Marrow","Lymph Node","Spleen","Liver","Lung", "Jejunum"),]

NK_prop$Tissue <- as.character(NK_prop$Tissue)
# Add number of data points to the label
for (my_tissue in unique(NK_prop$Tissue)){
  NK_prop$Tissue[NK_prop$Tissue == my_tissue] <- paste0(NK_prop$Tissue[NK_prop$Tissue == my_tissue], "\nn = ",
                                                        length(NK_prop$Tissue[NK_prop$Tissue == my_tissue]))
}
NK_prop$Tissue <- factor(NK_prop$Tissue, levels = c(NK_prop$Tissue[grep("Blood", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Bone Marrow", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Lymph Node", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Spleen", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Liver", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Lung", NK_prop$Tissue)[1]],
                                                    NK_prop$Tissue[grep("Jejunum", NK_prop$Tissue)[1]]))

# Make box plot
p <- ggplot(NK_prop, aes(x = Tissue, y = Percentage)) +
  geom_boxplot(fill = "firebrick", color = "black", size = 0.5, outlier.shape = NA, width = 0.6) + 
  stat_boxplot(geom = 'errorbar', size = 0.5, width = 0.3) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
  theme_classic() + labs(y = "NK Cell Percentage of Lymphocytes") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))

# Print to PDF
pdf(file = file.path(results_dir, "NK_proportion_of_lymphocytes.pdf"), width = 8, height = 4)
plot(p)
dev.off()



## Box plot of NK subset percentages across tissues ###################################################
# Load subset proportions data
NK_subset_proportions <- read.csv(file = file.path(data_dir, "NK_subset_proportions.csv"))

# Set AxLN tissues as LN for SIVS SAMPLES
NK_subset_proportions$Tissue[grepl("Ax_LN",NK_subset_proportions$Tissue)] <- "Lymph Node"

# Order factors
NK_subset_proportions$Tissue <- factor(NK_subset_proportions$Tissue, levels = c("Blood","Bone Marrow","Lymph Node","Spleen",
                                                                                "Liver","Lung", "Jejunum",
                                                                                "BAL", "Colon", "Ing_LN", "Mes_LN"))

# Only keep the tissues we want to show
NK_subset_proportions <- NK_subset_proportions[NK_subset_proportions$Tissue %in% c("Blood","Bone Marrow","Lymph Node","Spleen","Liver","Lung", "Jejunum"),]

NK_subset_proportions$Tissue <- as.character(NK_subset_proportions$Tissue)
# Add number of data points to the label
for (my_tissue in unique(NK_subset_proportions$Tissue)){
  NK_subset_proportions$Tissue[NK_subset_proportions$Tissue == my_tissue] <- paste0(NK_subset_proportions$Tissue[NK_subset_proportions$Tissue == my_tissue], "\nn = ",
                                                                                    length(NK_subset_proportions$Tissue[NK_subset_proportions$Tissue == my_tissue]))
}
NK_subset_proportions$Tissue <- factor(NK_subset_proportions$Tissue, levels = c(NK_subset_proportions$Tissue[grep("Blood", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Bone Marrow", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Lymph Node", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Spleen", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Liver", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Lung", NK_subset_proportions$Tissue)[1]],
                                                                                NK_subset_proportions$Tissue[grep("Jejunum", NK_subset_proportions$Tissue)[1]]))
# Melt into long formast
NK_subset_long <- reshape2::melt(NK_subset_proportions, id.vars = "Tissue", measure.vars = c("CD16_NK", "CD56_NK", "DN_NK"))

# Rename subsets
NK_subset_long$variable <- gsub("_", " ", NK_subset_long$variable)

# Order factors
NK_subset_long$variable <- gsub("CD56 NK", "CD56+ NK", NK_subset_long$variable)
NK_subset_long$variable <- gsub("CD16 NK", "CD16+ NK", NK_subset_long$variable)
NK_subset_long$variable <- factor(NK_subset_long$variable, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Make grouped box plot
p <- ggplot(NK_subset_long, aes(x = Tissue, y = value, fill = variable)) +
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3) +
  geom_point(position=position_jitterdodge(0.1), size = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  theme_classic() + labs(y = "Percentage of NK Cells") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"))

# Print to PDF
pdf(file = file.path(results_dir, "NK_subset_percentages.pdf"), width = 8, height = 4.6)
plot(p)
dev.off()



## Sankey plot of animal G45T peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'G45T_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


## Sankey plot of animal DGA0 peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'DGA0_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


# Repeat for DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# Repeat for CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


## Sankey plot of animal DGF5 peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'DGF5_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


# DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


## Sankey plot of animal TID peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'TID_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


## Bar plot of long term circulating proportions ###################################################
# Read the data
lt_circ_prop <- read.csv(file = file.path(data_dir, 'PB_long_term_circulating_proportions_with_T.csv'))

# Ensure the order of factors
lt_circ_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", lt_circ_prop$Cell_Type)
lt_circ_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", lt_circ_prop$Cell_Type)
lt_circ_prop$Cell_Type <- factor(lt_circ_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate standard error
lt_circ_prop_summary <- data_summary_se(lt_circ_prop, varname="percentage", 
                                        groupnames=c("Animal", "Cell_Type"))

# Ensure the order of factors
lt_circ_prop_summary$Animal <- factor(lt_circ_prop_summary$Animal, levels = c("DGA0 6hr", "DGF5 6hr","G45T 6hr", "TID 6hr", "G45T 48hr", "TID 48hr"))

# Replace space with new line
lt_circ_prop_summary$Animal <- gsub(" ", "\n-", lt_circ_prop_summary$Animal)
lt_circ_prop_summary$Animal <- factor(lt_circ_prop_summary$Animal, levels = c( "DGA0\n-6hr", "DGF5\n-6hr","G45T\n-6hr", "TID\n-6hr", "G45T\n-48hr", "TID\n-48hr"))

# Rename columns
colnames(lt_circ_prop_summary)[2] <- "Cell Type"
colnames(lt_circ_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(lt_circ_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(title = "", y = "Percentage Continuously Circulating") + 
  # scale_x_discrete(labels = levels(lt_circ_prop_summary$Animal)) +
  theme(axis.title.x  = element_blank(), axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"), legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "PB_long_term_circulating_percentages_with_T_cells.pdf"), width = 7*6/4, height = 4)
plot(p)
dev.off()


## Make dot plot comparing mean lt-circulating proportion of the 3 biological replicates ###################################################
# Load means from 6 hour positive proportions [means from previous plot normalized to lymphocyte percentage]
means_6hr <- read.csv(file = file.path(data_dir, "PB_6hr_means.csv"))

# Ensure the order of factors
means_6hr$Cell_Type <- gsub("CD56 NK", "CD56+ NK", means_6hr$Cell_Type)
means_6hr$Cell_Type <- gsub("CD16 NK", "CD16+ NK", means_6hr$Cell_Type)
means_6hr$Cell_Type <- factor(means_6hr$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Box plot with wilcoxon test to compare means
my_comparisons <- list( c("CD56+ NK", "DN NK"), c("DN NK", "CD16+ NK"), c("CD56+ NK", "CD16+ NK") )
p <- ggplot(means_6hr, aes(x = Cell_Type, y = value, fill = Cell_Type)) +
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F,
                     label.y = c(100, 110, 120)) +
  geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  labs(title = "", y = "Percentage Continuously Circulating") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")+
  scale_y_continuous(breaks=seq(0,100,25))


pdf(file = file.path(results_dir, "PB_6hour_proportions_boxplot.pdf"), width = 4, height = 4)
plot(p)
dev.off()


## Peripheral blood summary plot across animals ###################################################
# Read PB data
PB_data <- read.csv(file = file.path(data_dir, 'Summary_data_by_tissue_PB_with_T.csv'))

# Melt into long format
PB_data_long <- reshape2::melt(PB_data, id = c("Animal","Cell_Type"))

# Ensure the order of my factors
PB_data_long$Animal <- factor(PB_data_long$Animal, levels = c("DGA0", "DGF5", "G45T", "TID"))
PB_data_long$Cell_Type <- gsub("CD56 NK", "CD56+ NK", PB_data_long$Cell_Type)
PB_data_long$Cell_Type <- gsub("CD16 NK", "CD16+ NK", PB_data_long$Cell_Type)
PB_data_long$Cell_Type <- factor(PB_data_long$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16+ NK"))

# Rename Cell_Type column
colnames(PB_data_long)[2] <- "Cell Type"

# Rename the SIVS variable and ensure order
PB_data_long$variable <- gsub("SIVS_","", as.character(PB_data_long$variable))
PB_data_long$variable <- factor(PB_data_long$variable, levels = c("48hr", "24hr", "6hr", "2hr", "5min"))
PB_data_long$variable <- paste0("-",PB_data_long$variable)
PB_data_long$variable <- factor(PB_data_long$variable, levels = c("-48hr", "-24hr", "-6hr", "-2hr", "-5min"))

PB_data_long <- na.omit(PB_data_long)

# Make bar plot
p <- ggbarplot(PB_data_long, x = "variable", y = "value",
               add = c("mean_se", "jitter"),
               color = "Cell Type", shape = "Animal", fill = "grey95",
               size = 0.8,
               palette = brewer.pal(4, name = "Set1")[c(1,4,3,2)],
               position = position_dodge(0.8)) +
  labs(x = "SIVS Infusion", y = "Percentage SIVS Positive") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black")) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 0)
  )

# Save plot to PDF
pdf(file = file.path(results_dir, "PB_summary_figure_with_T.pdf"), width = 12, height = 5)
plot(p)
dev.off()


## Bar plot of IVAS+ percentages in lymph node ###################################################
# Read the data
LN_ivas_prop <- read.csv(file = file.path(data_dir, 'LN_IVAS_positive_proportions.csv'))

# # Exclude lymphocytes
# LN_ivas_prop <- LN_ivas_prop[LN_ivas_prop$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]

# Ensure the order of factors
LN_ivas_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_ivas_prop$Cell_Type)
LN_ivas_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_ivas_prop$Cell_Type)
LN_ivas_prop$Cell_Type <- factor(LN_ivas_prop$Cell_Type, levels = c("T Cell","CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate stadard deviation
LN_ivas_prop_summary <- data_summary_se(LN_ivas_prop, varname="ivas_percentage", 
                                        groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(LN_ivas_prop_summary)[2] <- "Cell Type"
colnames(LN_ivas_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_ivas_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "Percentage IVas+") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "LN_IVAS_positive_percentages_with_T.pdf"), width = 10*6/5, height = 4)
plot(p)
dev.off()

## LN IVAS+ summary dot plot 6hr ###################################################
# Load means normalized to T cell
means_ivas <- read.csv(file = file.path(data_dir, "LN_IVAS_positive_proportions_mean_dotplot.csv"))

# Ensure the order of factors
means_ivas$Cell_Type <- gsub("CD56 NK", "CD56+ NK", means_ivas$Cell_Type)
means_ivas$Cell_Type <- gsub("CD16 NK", "CD16+ NK", means_ivas$Cell_Type)
means_ivas$Cell_Type <- factor(means_ivas$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Box plot with wilcoxon test to compare means
my_comparisons <- list( c("CD56+ NK", "DN NK"), c("DN NK", "CD16+ NK"), c("CD56+ NK", "CD16+ NK") )
p <- ggplot(means_ivas, aes(x = Cell_Type, y = value, fill = Cell_Type)) +
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F,
                     label.y = c(log10(100), log10(150), log10(250))) +
  geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  scale_y_continuous(trans = "log10") +
  labs(title = "", y = "Percentage IVas+ (log scale)") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")

pdf(file = file.path(results_dir, "LN_IVAS_positive_boxplot.pdf"), width = 4, height = 4)
plot(p)
dev.off()


## Bar plot of long term tissue resident in lymph node ###################################################
# Read the data
LN_tr_prop <- read.csv(file = file.path(data_dir, 'LN_tissue_resident_proportions.csv'))

# # Exclude lymphocytes
# LN_tr_prop <- LN_tr_prop[LN_tr_prop$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]

# Ensure the order of factors
LN_tr_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_tr_prop$Cell_Type)
LN_tr_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_tr_prop$Cell_Type)
LN_tr_prop$Cell_Type <- factor(LN_tr_prop$Cell_Type, levels = c("T Cell","CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate stadard deviation
LN_tr_prop_summary <- data_summary_se(LN_tr_prop, varname="tr_percentage", 
                                      groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(LN_tr_prop_summary)[2] <- "Cell Type"
colnames(LN_tr_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_tr_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "Percentage Tissue Localized") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))


pdf(file = file.path(results_dir, "LN_tissue_resident_percentages_with_T.pdf"), width = 10*6/5, height = 4)
plot(p)
dev.off()


## LN tissue residence summary dot plot 6hr ###################################################
# Load means normalized to T cell
means_tr <- read.csv(file = file.path(data_dir, "LN_tissue_resident_proportions_mean_dotplot.csv"))

# Ensure the order of factors
means_tr$Cell_Type <- gsub("CD56 NK", "CD56+ NK", means_tr$Cell_Type)
means_tr$Cell_Type <- gsub("CD16 NK", "CD16+ NK", means_tr$Cell_Type)
means_tr$Cell_Type <- factor(means_tr$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Box plot with wilcoxon test to compare means
my_comparisons <- list( c("CD56+ NK", "DN NK"), c("DN NK", "CD16+ NK"), c("CD56+ NK", "CD16+ NK") )
p <- ggplot(means_tr, aes(x = Cell_Type, y = value, fill = Cell_Type)) +
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F,
                     label.y = c(105, 112.5, 120)) +
  geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  labs(title = "", y = "Percentage Tissue Localized") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,125), breaks=seq(0,100,25))

pdf(file = file.path(results_dir, "LN_tissue_residence_boxplot.pdf"), width = 4, height = 4)
plot(p)
dev.off()


## LN Entry rates bar plot 6hr ###################################################
# Read the data
LN_entry_prop <- read.csv(file = file.path(data_dir, 'LN_entry_rate_proportions_6hr.csv'))

# Ensure the order of factors
LN_entry_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_entry_prop$Cell_Type)
LN_entry_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_entry_prop$Cell_Type)
LN_entry_prop$Cell_Type <- factor(LN_entry_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate stadard deviation
LN_entry_prop_summary <- data_summary_se(LN_entry_prop, varname="entry_rate_6hr", 
                                         groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(LN_entry_prop_summary)[2] <- "Cell Type"
colnames(LN_entry_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_entry_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "Entry Rate into LN (% per Hour)") + 
  theme(axis.title.x  = element_blank(), axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"), legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "LN_entry_rates_6hr_with_T.pdf"), width = 10*6/5, height = 4)
plot(p)
dev.off()


## LN Entry rates summary dot plot 6hr ###################################################
# Load means from 6 hour LN entry rates [means from previous plot normalized to lymphocyte percentage]
means_LN_entry <- read.csv(file = file.path(data_dir, "LN_entry_rate_means_dotplot.csv"))

# Ensure the order of factors
means_LN_entry$Cell_Type <- gsub("CD56 NK", "CD56+ NK", means_LN_entry$Cell_Type)
means_LN_entry$Cell_Type <- gsub("CD16 NK", "CD16+ NK", means_LN_entry$Cell_Type)
means_LN_entry$Cell_Type <- factor(means_LN_entry$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Box plot with wilcoxon test to compare means
my_comparisons <- list( c("CD56+ NK", "DN NK"), c("DN NK", "CD16+ NK"), c("CD56+ NK", "CD16+ NK") )
p <- ggplot(means_LN_entry, aes(x = Cell_Type, y = value, fill = Cell_Type)) +
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F,
                     label.y = c(3.5, 3.75, 4)) +
  geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  labs(title = "", y = "Entry Rate into LN (% per Hour)") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")

pdf(file = file.path(results_dir, "LN_entry_rate_normalized_boxplot.pdf"), width = 4, height = 4)
plot(p)
dev.off()


## Lymph node summary plot across animals ###################################################
# Read PB data
LN_data <- read.csv(file = file.path(data_dir, 'Summary_data_by_tissue_LN_with_T.csv'))

# Melt into long format
LN_data_long <- reshape2::melt(LN_data, id = c("Animal_Sample","Cell_Type"))

# Ensure the order of my factors
LN_data_long$Animal <- factor(LN_data_long$Animal_Sample, levels = c("DGA0 AxLN", "DGA0 MesLN", "DGF5 AxLN", "G45T AxLN", "G45T IngLN", "TID AxLN"))
LN_data_long$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_data_long$Cell_Type)
LN_data_long$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_data_long$Cell_Type)
LN_data_long$Cell_Type <- factor(LN_data_long$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16+ NK"))

# Rename Cell_Type column
colnames(LN_data_long)[2] <- "Cell Type"

# Rename the SIVS variable and ensure order
LN_data_long$variable <- gsub("SIVS_","", as.character(LN_data_long$variable))
LN_data_long$variable <- factor(LN_data_long$variable, levels = c("48hr", "24hr", "6hr", "2hr", "5min"))
LN_data_long$variable <- paste0("-",LN_data_long$variable)
LN_data_long$variable <- factor(LN_data_long$variable, levels = c("-48hr", "-24hr", "-6hr", "-2hr", "-5min"))

LN_data_long <- na.omit(LN_data_long)

# Make bar plot
p <- ggbarplot(LN_data_long, x = "variable", y = "value",
               add = c("mean_se", "jitter"),
               color = "Cell Type", shape = "Animal", fill = "grey95",
               size = 0.8,
               palette = brewer.pal(4, name = "Set1")[c(1,4,3,2)],
               position = position_dodge(0.8)) +
  labs(x = "SIVS Infusion", y = "Percentage SIVS Positive") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black")) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 0)
  )

# Save plot to PDF
pdf(file = file.path(results_dir, "LN_summary_figure_with_T.pdf"), width = 12, height = 5)
plot(p)
dev.off()



## LN Phenotype plot for each animal ###################################################
# Read data
LN_pheno_data <- read.csv(file = file.path(data_dir, "LN_phenotype_data.csv"))

# Loop through each animal
for (my_animal in unique(LN_pheno_data$Animal)){
  LN_pheno_my_animal <- LN_pheno_data[LN_pheno_data$Animal == my_animal,]
  
  # Melt the DF
  LN_pheno_melt <- melt(LN_pheno_my_animal, id = c("Animal","Tissue"))
  
  # Remove the Ax_LN data
  LN_pheno_melt <- LN_pheno_melt[LN_pheno_melt$Tissue %in% c("PB", "IVAS_pos", "IVAS_neg"),]
  
  # Ensure order of factors
  LN_pheno_melt$Tissue <- factor(LN_pheno_melt$Tissue, levels = c("PB", "IVAS_pos", "IVAS_neg"))
  LN_pheno_melt$variable <- factor(LN_pheno_melt$variable, levels = c("CD56_NK", "DN_NK", "CD16_NK"))
  
  # Make phenotype pie chart for each subset
  pheno_plot_list <- list()
  mycols <- RColorBrewer::brewer.pal(4, "Set1")
  for (i in 1:length(unique(LN_pheno_melt$variable))){
    
    my_tissue <- unique(LN_pheno_melt$Tissue)[i]
    
    plot.data <- LN_pheno_melt[LN_pheno_melt$Tissue == my_tissue,]
    # Make pie chart
    pheno_plot_list[[i]] <- ggplot(plot.data,
                                   aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols[c(4,3,2)]) +
      theme_void() + theme(legend.position = "none")
    
  }
  
  # Make a PDF of the 3 pie charts for each tissue
  comb_plot_pheno <- cowplot::plot_grid(plotlist = pheno_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "NK_phenotype_pie_chart.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot_pheno)
  dev.off()
  
}


## LN Phenotype summary plot
# Read data
LN_pheno_data <- read.csv(file = file.path(data_dir, "LN_phenotype_data.csv"))

# Rename tissue variable
LN_pheno_data$Tissue <- plyr::revalue(LN_pheno_data$Tissue, c("PB" = "Blood", "IVAS_pos" = "LN IVas+", "IVAS_neg" = "LN IVas-") )

# Turn data into percentage out of 100
LN_pheno_data[,3:5] <- LN_pheno_data[,3:5] / rowSums(LN_pheno_data[,3:5]) * 100

# Melt
LN_pheno_melt <- reshape2::melt(LN_pheno_data, id = c("Animal","Tissue"))

# Replace _ with space in subset name
LN_pheno_melt$variable <- gsub("_", " ", LN_pheno_melt$variable)
LN_pheno_melt$variable <- gsub("CD56 NK", "CD56+ NK", LN_pheno_melt$variable)
LN_pheno_melt$variable <- gsub("CD16 NK", "CD16+ NK", LN_pheno_melt$variable)
LN_pheno_melt$variable <- factor(LN_pheno_melt$variable, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

# Keep only the desired data
LN_pheno_melt <- LN_pheno_melt[LN_pheno_melt$Tissue %in% c("Blood", "LN IVas+", "LN IVas-"),]

# Ensure order of factors
LN_pheno_melt$Tissue <- factor(LN_pheno_melt$Tissue, levels = c("Blood", "LN IVas+", "LN IVas-"))

# Set comparisos list for statistics
my_comparisons <- list( c("Blood", "LN IVas+"), c("LN IVas+", "LN IVas-"), c("Blood", "LN IVas-") )

# Box plot with wilcoxon test to compare means
p <- ggplot(LN_pheno_melt, aes(x = Tissue, y = value)) +
  facet_wrap(~variable)+
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6, aes(fill = variable)) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F,
                     label.y = c(100, 107, 114)) +
  # geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  geom_jitter(size = 2.5, position=position_dodge(0.8)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  scale_y_continuous(trans = "log10") +
  labs(title = "", y = "Percentage of NK cells") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12,  color = "black", angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 12, color = "black", face = "bold")) +
  scale_y_continuous(breaks=seq(0,100,20))

# Save to PDF
pdf(file = file.path(results_dir, "LN_phenotype_summary_boxplot.pdf"), width = 7, height = 5)
plot(p)
dev.off()


## Spleen Phenotype plot for each animal ###################################################
# Read data
tissue_pheno_data <- read.csv(file = file.path(data_dir, "Spleen_phenotype_data.csv"))

# Loop through each animal
for (my_animal in unique(tissue_pheno_data$Animal)){
  pheno_my_animal <- tissue_pheno_data[tissue_pheno_data$Animal == my_animal,]
  
  # Melt the DF
  pheno_melt <- melt(pheno_my_animal, id = c("Animal","Tissue"))
  
  # Select data
  pheno_melt <- pheno_melt[pheno_melt$Tissue %in% c("PB", "Spleen IVAS_pos", "Spleen IVAS_neg"),]
  
  # Ensure order of factors
  pheno_melt$Tissue <- factor(pheno_melt$Tissue, levels = c("PB", "Spleen IVAS_pos", "Spleen IVAS_neg"))
  pheno_melt$variable <- factor(pheno_melt$variable, levels = c("CD56_NK", "DN_NK", "CD16_NK"))
  
  # Make phenotype pie chart for each subset
  pheno_plot_list <- list()
  mycols <- RColorBrewer::brewer.pal(4, "Set1")
  for (i in 1:length(unique(pheno_melt$variable))){
    
    my_tissue <- unique(pheno_melt$Tissue)[i]
    
    plot.data <- pheno_melt[pheno_melt$Tissue == my_tissue,]
    # Make pie chart
    pheno_plot_list[[i]] <- ggplot(plot.data,
                                   aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols[c(4,3,2)]) +
      theme_void() + theme(legend.position = "none")
    
  }
  
  # Make a PDF of the 3 pie charts for each tissue
  comb_plot_pheno <- cowplot::plot_grid(plotlist = pheno_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Spleen_NK_phenotype_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot_pheno)
  dev.off()
  
}


## Bar plot of IVAS+ percentages in Spleen ###################################################
# Read the data
Spleen_ivas_prop <- read.csv(file = file.path(data_dir, 'Spleen_IVAS_positive_proportions_with_T.csv'))

# Exclude lymphocytes
Spleen_ivas_prop <- Spleen_ivas_prop[Spleen_ivas_prop$Cell_Type %in% c("T Cell", "CD56 NK", "DN NK", "CD16 NK"),]

# Ensure the order of factors
Spleen_ivas_prop$Cell_Type <- factor(Spleen_ivas_prop$Cell_Type, levels = c("T Cell","CD56 NK", "DN NK", "CD16 NK"))

# Calculate stadard deviation
Spleen_ivas_prop_summary <- data_summary_se(Spleen_ivas_prop, varname="ivas_percentage", 
                                            groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(Spleen_ivas_prop_summary)[2] <- "Cell Type"
colnames(Spleen_ivas_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(Spleen_ivas_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "Percentage IVas+") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "Spleen_IVAS_positive_percentages_with_T.pdf"), width = 6, height = 4)
plot(p)
dev.off()



## Bar plot of long term tissue resident in Spleen ###################################################
# Read the data
Spleen_tr_prop <- read.csv(file = file.path(data_dir, 'Spleen_tissue_resident_proportions_with_T.csv'))

# Exclude lymphocytes
Spleen_tr_prop <- Spleen_tr_prop[Spleen_tr_prop$Cell_Type %in% c("T Cell", "CD56 NK", "DN NK", "CD16 NK"),]

# Ensure the order of factors
Spleen_tr_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", Spleen_tr_prop$Cell_Type)
Spleen_tr_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", Spleen_tr_prop$Cell_Type)
Spleen_tr_prop$Cell_Type <- factor(Spleen_tr_prop$Cell_Type, levels = c("T Cell","CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate stadard deviation
Spleen_tr_prop_summary <- data_summary_se(Spleen_tr_prop, varname="tr_percentage", 
                                          groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(Spleen_tr_prop_summary)[2] <- "Cell Type"
colnames(Spleen_tr_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(Spleen_tr_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "Percentage Tissue Localized") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "Spleen_tissue_resident_percentages_with_T.pdf"), width = 6, height = 4)
plot(p)
dev.off()



## Spleen Entry rates bar plot 6hr ###################################################
# Read the data
Spleen_entry_prop <- read.csv(file = file.path(data_dir, 'Spleen_entry_rate_proportions_6hr_updated.csv'))

# Ensure the order of factors
Spleen_entry_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", Spleen_entry_prop$Cell_Type)
Spleen_entry_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", Spleen_entry_prop$Cell_Type)
Spleen_entry_prop$Cell_Type <- factor(Spleen_entry_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16+ NK"))

# Calculate stadard deviation
Spleen_entry_prop_summary <- data_summary_se(Spleen_entry_prop, varname="entry_rate_6hr", 
                                             groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(Spleen_entry_prop_summary)[2] <- "Cell Type"
colnames(Spleen_entry_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(Spleen_entry_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3,2)]) + labs(y = "IVAS- Entry Rate into Spleen (% per Hour)") + 
  theme(axis.title.x  = element_blank(), axis.title.y  = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"), legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "Spleen_entry_rates_6hr_with_T.pdf"), width = 6, height = 4)
plot(p)
dev.off()



## Liver Composition plot for each animal ###################################################
# Read data
tissue_comp_data <- read.csv(file = file.path(data_dir, "Liver_composition_data.csv"))

# Loop through each animal
for (my_animal in unique(tissue_comp_data$Animal)){
  tissue_comp_my_animal <- tissue_comp_data[tissue_comp_data$Animal == my_animal,]
  
  # Melt the DF
  comp_my_tissue_melt <- melt(tissue_comp_my_animal, id = c("Animal", "Tissue", "Cell_Type"))
  
  # Exclude lymphocytes
  comp_my_tissue_melt <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]
  
  # Ensure order of factors
  comp_my_tissue_melt$Cell_Type <- factor(comp_my_tissue_melt$Cell_Type, levels = c("CD56 NK", "DN NK", "CD16 NK"))
  
  # Make IVAS+/- pie chart for each subset
  IVAS_plot_list <- list()
  TR_plot_list <- list()
  mycols1 <- RColorBrewer::brewer.pal(3, "Set1")
  mycols2 <- RColorBrewer::brewer.pal(9, "Blues")
  for (i in 1:length(levels(comp_my_tissue_melt$Cell_Type))){
    
    my_subset <- levels(comp_my_tissue_melt$Cell_Type)[i]
    
    # Make IVAS+/- pie chart
    plot.data <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type == my_subset &
                                       comp_my_tissue_melt$variable %in% c("IVAS_pos", "IVAS_neg"),]
    
    IVAS_plot_list[[i]] <- ggplot(plot.data,
                                  aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols1[1:2]) +
      theme_void() + theme(legend.position = "none")
    
    
  }
  
  # Make a PDF of the IVAS+/- pie charts for each tissue
  comb_plot1 <- cowplot::plot_grid(plotlist = IVAS_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Liver_IVAS_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot1)
  dev.off()
  
}


## Lung Composition plot for each animal ###################################################
# Read data
tissue_comp_data <- read.csv(file = file.path(data_dir, "Lung_composition_data.csv"))

# Loop through each animal
for (my_animal in unique(tissue_comp_data$Animal)){
  tissue_comp_my_animal <- tissue_comp_data[tissue_comp_data$Animal == my_animal,]
  
  # Melt the DF
  comp_my_tissue_melt <- melt(tissue_comp_my_animal, id = c("Animal", "Tissue", "Cell_Type"))
  
  # Exclude lymphocytes
  comp_my_tissue_melt <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]
  
  # Ensure order of factors
  comp_my_tissue_melt$Cell_Type <- factor(comp_my_tissue_melt$Cell_Type, levels = c("CD56 NK", "DN NK", "CD16 NK"))
  
  # Make IVAS+/- pie chart for each subset
  IVAS_plot_list <- list()
  TR_plot_list <- list()
  mycols1 <- RColorBrewer::brewer.pal(3, "Set1")
  mycols2 <- RColorBrewer::brewer.pal(9, "Blues")
  for (i in 1:length(levels(comp_my_tissue_melt$Cell_Type))){
    
    my_subset <- levels(comp_my_tissue_melt$Cell_Type)[i]
    
    # Make IVAS+/- pie chart
    plot.data <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type == my_subset &
                                       comp_my_tissue_melt$variable %in% c("IVAS_pos", "IVAS_neg"),]
    
    IVAS_plot_list[[i]] <- ggplot(plot.data,
                                  aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols1[1:2]) +
      theme_void() + theme(legend.position = "none")
    
    # Make Tissue residence pie chart
    plot.data <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type == my_subset &
                                       comp_my_tissue_melt$variable %in% c("RI", "TR"),]
    
    TR_plot_list[[i]] <- ggplot(plot.data,
                                aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols2[c(5,9)]) +
      theme_void() + theme(legend.position = "none")
    
  }
  
  # Make a PDF of the IVAS+/- pie charts for each tissue
  comb_plot1 <- cowplot::plot_grid(plotlist = IVAS_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Lung_IVAS_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot1)
  dev.off()
  
  
  # Make a PDF of the TR pie charts for each tissue 
  comb_plot2 <- cowplot::plot_grid(plotlist = TR_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Lung_TR_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot2)
  dev.off()
  
}


## Jejunum Composition plot for each animal ###################################################
# Read data
tissue_comp_data <- read.csv(file = file.path(data_dir, "Jejunum_composition_data.csv"))

# Loop through each animal
for (my_animal in unique(tissue_comp_data$Animal)){
  tissue_comp_my_animal <- tissue_comp_data[tissue_comp_data$Animal == my_animal,]
  
  # Melt the DF
  comp_my_tissue_melt <- melt(tissue_comp_my_animal, id = c("Animal", "Tissue", "Cell_Type"))
  
  # Exclude lymphocytes
  comp_my_tissue_melt <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]
  
  # Ensure order of factors
  comp_my_tissue_melt$Cell_Type <- factor(comp_my_tissue_melt$Cell_Type, levels = c("CD56 NK", "DN NK", "CD16 NK"))
  
  # Make IVAS+/- pie chart for each subset
  IVAS_plot_list <- list()
  TR_plot_list <- list()
  mycols1 <- RColorBrewer::brewer.pal(3, "Set1")
  mycols2 <- RColorBrewer::brewer.pal(9, "Blues")
  for (i in 1:length(levels(comp_my_tissue_melt$Cell_Type))){
    
    my_subset <- levels(comp_my_tissue_melt$Cell_Type)[i]
    
    # Make IVAS+/- pie chart
    plot.data <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type == my_subset &
                                       comp_my_tissue_melt$variable %in% c("IVAS_pos", "IVAS_neg"),]
    
    IVAS_plot_list[[i]] <- ggplot(plot.data,
                                  aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols1[1:2]) +
      theme_void() + theme(legend.position = "none")
    
    # Make Tissue residence pie chart
    plot.data <- comp_my_tissue_melt[comp_my_tissue_melt$Cell_Type == my_subset &
                                       comp_my_tissue_melt$variable %in% c("RI", "TR"),]
    
    TR_plot_list[[i]] <- ggplot(plot.data,
                                aes(x = "", y = value, fill = variable)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #   geom_text(aes(y = lab.ypos, label = variable), color = "black")+
      scale_fill_manual(values = mycols2[c(5,9)]) +
      theme_void() + theme(legend.position = "none")
    
  }
  
  # Make a PDF of the IVAS+/- pie charts for each tissue
  comb_plot1 <- cowplot::plot_grid(plotlist = IVAS_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Jejunum_IVAS_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot1)
  dev.off()
  
  
  # Make a PDF of the TR pie charts for each tissue 
  comb_plot2 <- cowplot::plot_grid(plotlist = TR_plot_list, align = "hv", ncol = 3)
  pdf(file = file.path(results_dir, paste(my_animal, "Jejunum_TR_pie_chart_updated.pdf", sep = "_")), width = 5*3/2, height = 2.5)
  plot(comb_plot2)
  dev.off()
  
}



## Comparison of each surface marker between LN RI and TR
# Read data
tissue_comp_data_full <- read.csv(file = file.path(data_dir, "RI_vs_TR_LN_marker_percentages.csv"))

# Loop through each marker
for (i in 1:length(unique(tissue_comp_data_full$Marker))){
  my_marker <- unique(tissue_comp_data_full$Marker)[i]
  tissue_comp_data <- tissue_comp_data_full[tissue_comp_data_full$Marker == my_marker,]
  
  # Set order
  tissue_comp_data$Cell_Type <- gsub("CD56 NK", "CD56+ NK", tissue_comp_data$Cell_Type)
  tissue_comp_data$Cell_Type <- gsub("CD16 NK", "CD16+ NK", tissue_comp_data$Cell_Type)
  tissue_comp_data$Cell_Type <- factor(tissue_comp_data$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))
  
  tissue_comp_data$Subset <- gsub("TR","TL", tissue_comp_data$Subset)
  
  # Set comparisos list for statistics
  my_comparisons <- list( c("RI", "TL"))
  

  # Box plot with wilcoxon test to compare means
  p <- ggplot(tissue_comp_data, aes(x = Subset, y = percentage)) +
    facet_wrap(~Cell_Type)+
    geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6, aes(fill = Cell_Type)) + 
    stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 10) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F) +
    # geom_point(position=position_jitterdodge(0.1), size = 2.5) +
    geom_jitter(size = 2.5, position=position_dodge(0.8)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
    labs(title = "", y = paste0("Percentage ", my_marker, "+")) + 
    theme_classic() + 
    theme(title = element_blank(),
          plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
          axis.title.y = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12,  color = "black", angle = 45, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 12, color = "black", face = "bold"))
  
  # Save to PDF
  pdf(file = file.path(results_dir, paste0("RI_vs_TR_",my_marker,"_percentage_boxplot.pdf")), width = 6, height = 5)
  plot(p)
  dev.off()
}



## Comparison of percentage TL between CD69+ and CD69- from LN
# Read data
tissue_comp_data <- read.csv(file = file.path(data_dir, "LN_CD69_pos_neg_TL_percentages.csv"))

# Set order
tissue_comp_data$Cell_Type <- gsub("CD56 NK", "CD56+ NK", tissue_comp_data$Cell_Type)
tissue_comp_data$Cell_Type <- gsub("CD16 NK", "CD16+ NK", tissue_comp_data$Cell_Type)
tissue_comp_data$Cell_Type <- factor(tissue_comp_data$Cell_Type, levels = c("CD56+ NK", "DN NK", "CD16+ NK"))

tissue_comp_data$Subset <- gsub("CD69_pos","CD69+", tissue_comp_data$Subset)
tissue_comp_data$Subset <- gsub("CD69_neg","CD69-", tissue_comp_data$Subset)

# Set comparisos list for statistics
my_comparisons <- list( c("CD69+", "CD69-"))

# Box plot with wilcoxon test to compare means
p <- ggplot(tissue_comp_data, aes(x = Subset, y = percentage)) +
  facet_wrap(~Cell_Type)+
  geom_boxplot(color = "black", size = 0.4, outlier.shape = NA, position=position_dodge(0.75), width = 0.6, aes(fill = Cell_Type)) + 
  stat_boxplot(position=position_dodge(0.75), geom = 'errorbar', size = 0.5, width = 0.3, coef = 10) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F, label.y = 102) +
  # geom_point(position=position_jitterdodge(0.1), size = 2.5) +
  geom_jitter(size = 2.5, position=position_dodge(0.8)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[c(4,3,2)]) +
  labs(title = "", y = "Percentage Tissue Localized") + 
  theme_classic() + 
  theme(title = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)),
        axis.title.y = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12,  color = "black", angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 12, color = "black", face = "bold")) +
  scale_y_continuous(breaks=seq(0,100,25), limits = c(0,107))

# Save to PDF
pdf(file = file.path(results_dir, paste0("CD69_pos_vs_neg_TL_percentages.pdf")), width = 6, height = 5)
plot(p)
dev.off()




## Sankey plot of animal HAXR peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'HAXR_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]
links$source <- as.factor(links$source)
links$target <- as.factor(links$target)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16dim NK
links <- sankey_data[sankey_data$cell_type == "CD16dim NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p



## Sankey plot of animal K628 peripheral blood data ###################################################
# Read data
sankey_data <- read.csv(file = file.path(data_dir,'K628_PB_data_for_Sankey.csv'))

# CD56 NK
links <- sankey_data[sankey_data$cell_type == "CD56 NK", 2:4]
links$source <- as.factor(links$source)
links$target <- as.factor(links$target)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Purples")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#54278F", "#9E9AC8"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# DN NK
links <- sankey_data[sankey_data$cell_type == "DN NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Greens")[c(8,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#006D2C", "#74C476"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16dim NK
links <- sankey_data[sankey_data$cell_type == "CD16dim NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p

# CD16 NK
links <- sankey_data[sankey_data$cell_type == "CD16 NK", 2:4]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Color nodes
nodes$group <- as.factor(c("a","b","a","b","a","b"))
brewer.pal(9, "Blues")[c(9,5)]
my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#08306B", "#6BAED6"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize = 18,
                   height = 300, width = 450,
                   colourScale = my_color, NodeGroup = "group")
p


 
## Bar plot of long term circulating proportions comparison between depletion and SS ########################################
# Read the data
lt_circ_prop <- read.csv(file = file.path(data_dir, 'PB_long_term_circulating_proportions_CD16_depletion.csv'))

# Ensure the order of factors
lt_circ_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", lt_circ_prop$Cell_Type)
lt_circ_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", lt_circ_prop$Cell_Type)
lt_circ_prop$Cell_Type <- factor(lt_circ_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16dim NK", "CD16+ NK"))

# Calculate standard error
lt_circ_prop_summary <- data_summary_se(lt_circ_prop, varname="percentage", 
                                        groupnames=c("Animal", "Cell_Type"))

# Ensure the order of factors
lt_circ_prop_summary$Animal <- factor(lt_circ_prop_summary$Animal, levels = c("G45T 24hr", "TID 24hr", "G45T 48hr", "TID 48hr", "HAXR 24hr",  "K628 24hr", "HAXR 48hr"))

# Replace space with new line
lt_circ_prop_summary$Animal <- gsub(" ", "\n-", lt_circ_prop_summary$Animal)
lt_circ_prop_summary$Animal <- factor(lt_circ_prop_summary$Animal, levels = c( "G45T\n-24hr", "TID\n-24hr", "G45T\n-48hr", "TID\n-48hr", "HAXR\n-24hr", "K628\n-24hr", "HAXR\n-48hr"))


# Rename columns
colnames(lt_circ_prop_summary)[2] <- "Cell Type"
colnames(lt_circ_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(lt_circ_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3)], RColorBrewer::brewer.pal(11, "RdBu")[c(8)], RColorBrewer::brewer.pal(4, "Set1")[c(2)])) +
  labs(title = "", y = "Percentage Continuously Circulating") + 
  # scale_x_discrete(labels = levels(lt_circ_prop_summary$Animal)) +
  theme(axis.title.x  = element_blank(), axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"), legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "PB_long_term_circulating_percentages_CD16_depletion.pdf"), width = 8.5*7/4, height = 4)
plot(p)
dev.off()


## Bar plot of IVAS+ percentages in lymph node comparing HAXR and SS mean ###################################################
# Read the data
LN_ivas_prop <- read.csv(file = file.path(data_dir, 'LN_IVAS_positive_proportions_CD16_depletion.csv'))

# # Exclude lymphocytes
# LN_ivas_prop <- LN_ivas_prop[LN_ivas_prop$Cell_Type %in% c("CD56 NK", "DN NK", "CD16 NK"),]

# Ensure the order of factors
LN_ivas_prop$Group <- factor(LN_ivas_prop$Group, levels = c("AxLN Mean", "HAXR AxLN", "K628 AxLN"))
LN_ivas_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_ivas_prop$Cell_Type)
LN_ivas_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_ivas_prop$Cell_Type)
LN_ivas_prop$Cell_Type <- factor(LN_ivas_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16dim NK", "CD16+ NK"))
# LN_ivas_prop$Cell_Type <- factor(LN_ivas_prop$Cell_Type, levels = c("T Cell","CD56+ NK", "DN NK", "CD16 Dim", "CD16 Bright"))

# Calculate stadard deviation
LN_ivas_prop_summary <- data_summary_se(LN_ivas_prop, varname="ivas_percentage", 
                                        groupnames=c("Group", "Cell_Type"))

# Rename colums
colnames(LN_ivas_prop_summary)[2] <- "Cell Type"
colnames(LN_ivas_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_ivas_prop_summary, aes(x=Group, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3)], RColorBrewer::brewer.pal(11, "RdBu")[c(8)], RColorBrewer::brewer.pal(4, "Set1")[c(2)])) +
  labs(y = "Percentage IVas+") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "LN_IVAS_positive_percentages_CD16_depletion.pdf"), width = 8, height = 4)
plot(p)
dev.off()



## HAXR Bar plot of long term tissue resident in lymph node vs mean S.S. ###################################################
# Read the data
LN_tr_prop <- read.csv(file = file.path(data_dir, 'LN_tissue_resident_proportions_CD16_depletion.csv'))

# Include only CD56 and DN NK cells
# LN_tr_prop <- LN_tr_prop[LN_tr_prop$Cell_Type %in% c("T Cell", "CD56 NK", "DN NK"),]

# Ensure the order of factors
LN_tr_prop$Group <- factor(LN_tr_prop$Group, levels = c("AxLN Mean", "HAXR AxLN", "K628 AxLN"))
LN_tr_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_tr_prop$Cell_Type)
LN_tr_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_tr_prop$Cell_Type)
LN_tr_prop$Cell_Type <- factor(LN_tr_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16dim NK", "CD16+ NK"))
# LN_tr_prop$Cell_Type <- factor(LN_tr_prop$Cell_Type, levels = c("T Cell","CD56 NK", "DN NK", "CD16 Dim", "CD16 Bright"))

# Calculate stadard deviation
LN_tr_prop_summary <- data_summary_se(LN_tr_prop, varname="tr_percentage", 
                                      groupnames=c("Group", "Cell_Type"))

# Rename colums
colnames(LN_tr_prop_summary)[2] <- "Cell Type"
colnames(LN_tr_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_tr_prop_summary, aes(x=Group, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3)], RColorBrewer::brewer.pal(11, "RdBu")[c(8)], RColorBrewer::brewer.pal(4, "Set1")[c(2)])) +
  labs(y = "Percentage Tissue Localized") + 
  theme(axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "LN_tissue_localized_percentages_CD16_depletion.pdf"), width =8, height = 4)
plot(p)
dev.off()


## LN Entry rates bar plot 48hr ###################################################
# Read the data
LN_entry_prop <- read.csv(file = file.path(data_dir, 'LN_entry_rate_proportions_CD16_depletion.csv'))

# Include only CD56 and DN NK cells
# LN_entry_prop <- LN_entry_prop[LN_entry_prop$Cell_Type %in% c("T Cell", "CD56 NK", "DN NK"),]

# Ensure the order of factors
LN_entry_prop$Cell_Type <- gsub("CD56 NK", "CD56+ NK", LN_entry_prop$Cell_Type)
LN_entry_prop$Cell_Type <- gsub("CD16 NK", "CD16+ NK", LN_entry_prop$Cell_Type)
# LN_entry_prop$Cell_Type <- gsub("CD16dim NK", "CD16 Dim", LN_entry_prop$Cell_Type)
# LN_entry_prop$Cell_Type <- factor(LN_entry_prop$Cell_Type, levels = c("T Cell", "CD56+ NK", "DN NK", "CD16 Dim", "CD16 Bright"))
LN_entry_prop$Cell_Type <- factor(LN_entry_prop$Cell_Type, levels = c("T Cell","CD56+ NK", "DN NK", "CD16dim NK", "CD16+ NK"))

# Include only -48 hour calculated rates
LN_entry_prop <- LN_entry_prop[LN_entry_prop$Animal %in% c("G45T AxLN -24hr", "TID AxLN -24hr", "HAXR AxLN -24hr", "K628 AxLN -24hr"),]

# Add line break
LN_entry_prop$Animal <- gsub(" -", "\n-", LN_entry_prop$Animal)
LN_entry_prop$Animal <- factor(LN_entry_prop$Animal, levels = c( "G45T AxLN\n-24hr", "TID AxLN\n-24hr", "HAXR AxLN\n-24hr","K628 AxLN\n-24hr"))

# Calculate stadard deviation
LN_entry_prop_summary <- data_summary_se(LN_entry_prop, varname="entry_rate", 
                                         groupnames=c("Animal", "Cell_Type"))

# Rename colums
colnames(LN_entry_prop_summary)[2] <- "Cell Type"
colnames(LN_entry_prop_summary)[3] <- "value"

# Bar chart with error bars on top
p <- ggplot(LN_entry_prop_summary, aes(x=Animal, y=value, fill=`Cell Type`)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.75) +
  geom_errorbar(aes(ymin=value, ymax=value+se), width=.2,
                position=position_dodge(.75)) + theme_classic() +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Set1")[c(1,4,3)], RColorBrewer::brewer.pal(11, "RdBu")[c(8)], RColorBrewer::brewer.pal(4, "Set1")[c(2)])) +
  labs(y = "Entry Rate into LN (% per Hour)") + 
  theme(axis.title.x  = element_blank(), axis.title.y  = element_text(size = 14, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"), legend.text = element_text(size = 14, color = "black"))

pdf(file = file.path(results_dir, "LN_entry_rates_24hr_CD16_depletion.pdf"), width = 10, height = 4)
plot(p)
dev.off()


