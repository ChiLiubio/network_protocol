

# show other distance matrix results except for phylogenetic distance

if(!require(agricolae)){
	install.packages("agricolae")
}

# load microeco package
library(microeco)
# load magrittr package to use pipe operator
library(magrittr)
set.seed(123)
library(meconetcomp)

library(ggplot2)
theme_set(theme_bw())



data(stool_met)

load("stool_met_network.RData")


############################
# Niche overlap

# require MicroNiche
if(!require(MicroNiche)){
	install.packages("MicroNiche")
}
library(MicroNiche)

tmp <- stool_met$otu_table
tmp <- data.frame(Taxon = rownames(tmp), tmp)
# there is a bug in levins.overlap of MicroNiche 
tmp <- rbind(tmp[1, ], tmp)
tmp$Taxon %<>% as.character
tmp[1, 1] <- "add"
tmp$Taxon %<>% as.factor

levins_overlap <- levins.overlap(tmp)
rownames(levins_overlap) <- levins_overlap[, 1]
levins_overlap %<>% .[, -1] %>% as.matrix

tmp <- edge_node_distance$new(network_list = stool_met_network, dis_matrix = levins_overlap, label = c("+", "-"))
tmp$cal_diff(method = "anova")
tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5, xtext_angle = 30)


##########################################
# Competition index from genome analysis

load("genome_interactions.RData")

tmp <- edge_node_distance$new(network_list = stool_met_network, dis_matrix = genome_interactions$competition, label = c("+", "-"))
tmp$cal_diff(method = "anova")
tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5, xtext_angle = 30) + ylab("Competition index")


#########################
# Complementarity index

tmp <- edge_node_distance$new(network_list = stool_met_network, dis_matrix = genome_interactions$complementarity, label = c("+", "-"))
tmp$cal_diff(method = "anova")
tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5, xtext_angle = 30) + ylab("Complementarity index")






