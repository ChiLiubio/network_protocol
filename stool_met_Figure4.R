
# load microeco package
library(microeco)
# load magrittr package to use pipe operator
library(magrittr)
set.seed(123)
library(tibble)
library(dplyr)

library(ggplot2)
theme_set(theme_bw())


data("stool_met.RData")
data("stool_met_network.RData")
source("Funtion_utilities.R")

############################
# Niche overlap

# require MicroNiche
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

res <- data.frame()
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "+", select_matrix = levins_overlap, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "-", select_matrix = levins_overlap, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}

res$Group %<>% factor(., levels = names(stool_met_network))
res %<>% .[!is.na(.$Value), ]
res$Measure <- res$Group
res$Group <- paste0(res$Measure, "_", res$label)
# obtain significance
tmp1 <- res
tmp1$Measure <- "Test"
tmp2 <- trans_alpha$new(dataset = NULL)
tmp2$data_alpha <- tmp1
tmp2$group <- "Group"
tmp2$cal_diff(method = "anova", measure = "Test")

tmp3 <- tmp2$res_diff %>% tibble::rownames_to_column(var = "Group")
tmp5 <- res[, c(1, 3, 4)] %>% unique
tmp5 <- left_join(tmp3, tmp5, by = c("Group" = "Group"))
colnames(tmp5)[colnames(tmp5) == "Group"] <- "Type"
colnames(tmp5)[colnames(tmp5) == "label"] <- "Group"
colnames(tmp5)[colnames(tmp5) == "Test"] <- "Significance"

colnames(res)[colnames(res) == "Group"] <- "Type"
colnames(res)[colnames(res) == "label"] <- "Group"
res$Group %<>% factor(., levels = c("+", "-"))
res$Measure %<>% factor(., levels = names(stool_met_network))

g1 <- ggplot(res, aes(x = Measure, y = Value, color = Group)) +
	theme_bw() +
	geom_boxplot(position = position_dodge(0.9)) +
	geom_text(aes(x = x, y = y, label = add), data = get_sig_data(abund_data = res, diff_data = tmp5), inherit.aes = FALSE) +
	scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
	ylab("Niche overlap") +
	xlab("") +
	theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 12)) +
	theme(axis.title.y = element_text(size = 18)) +
	theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm")) + 
	labs(color = "Label")

g1


##########################################
# Phylogenetic distance

phylogenetic_distance <- cophenetic(stool_met$phylo_tree) %>% as.matrix

res <- data.frame()
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "+", select_matrix = phylogenetic_distance, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "-", select_matrix = phylogenetic_distance, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}

res$Group %<>% factor(., levels = names(stool_met_network))
res %<>% .[!is.na(.$Value), ]
res$Measure <- res$Group
res$Group <- paste0(res$Measure, "_", res$label)
# obtain significance
tmp1 <- res
tmp1$Measure <- "Test"
tmp2 <- trans_alpha$new(dataset = NULL)
tmp2$data_alpha <- tmp1
tmp2$group <- "Group"
tmp2$cal_diff(method = "anova", measure = "Test")

tmp3 <- tmp2$res_diff
tmp4 <- tmp3 %>% rownames_to_column(var = "Group")
tmp5 <- res[, c(1, 3, 4)] %>% unique
tmp5 <- left_join(tmp4, tmp5, by = c("Group" = "Group"))
colnames(tmp5)[colnames(tmp5) == "Group"] <- "Type"
colnames(tmp5)[colnames(tmp5) == "label"] <- "Group"
colnames(tmp5)[colnames(tmp5) == "Test"] <- "Significance"

colnames(res)[colnames(res) == "Group"] <- "Type"
colnames(res)[colnames(res) == "label"] <- "Group"
res$Group %<>% factor(., levels = c("+", "-"))
res$Measure %<>% factor(., levels = names(stool_met_network))


g2 <- ggplot(res, aes(x = Measure, y = Value, color = Group)) +
	theme_bw() +
	geom_boxplot(position = position_dodge(0.9)) +
	geom_text(aes(x = x, y = y, label = add), data = get_sig_data(abund_data = res, diff_data = tmp5), inherit.aes = FALSE) +
	scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
	ylab("Phylogenetic distance") +
	xlab("") +
	theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 12)) +
	theme(axis.title.y = element_text(size = 18)) +
	theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm")) + 
	labs(color = "Label")

g2


##########################################
# Competition index

load("genome_interactions.RData")

x1 <- genome_interactions$competition

res <- data.frame()
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "+", select_matrix = x1, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "-", select_matrix = x1, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}

res$Group %<>% factor(., levels = names(stool_met_network))
res %<>% .[!is.na(.$Value), ]
res$Measure <- res$Group
res$Group <- paste0(res$Measure, "_", res$label)
# Significance
tmp1 <- res
tmp1$Measure <- "Test"
tmp2 <- trans_alpha$new(dataset = NULL)
tmp2$data_alpha <- tmp1
tmp2$group <- "Group"
tmp2$cal_diff(method = "anova", measure = "Test")

tmp3 <- tmp2$res_diff
tmp4 <- tmp3 %>% rownames_to_column(var = "Group")
tmp5 <- res[, c(1, 3, 4)] %>% unique
tmp5 <- left_join(tmp4, tmp5, by = c("Group" = "Group"))
colnames(tmp5)[colnames(tmp5) == "Group"] <- "Type"
colnames(tmp5)[colnames(tmp5) == "label"] <- "Group"
colnames(tmp5)[colnames(tmp5) == "Test"] <- "Significance"

colnames(res)[colnames(res) == "Group"] <- "Type"
colnames(res)[colnames(res) == "label"] <- "Group"
res$Group %<>% factor(., levels = c("+", "-"))
res$Measure %<>% factor(., levels = names(stool_met_network))


g3 <- ggplot(res, aes(x = Measure, y = Value, color = Group)) +
	theme_bw() +
	geom_boxplot(position = position_dodge(0.9)) +
	geom_text(aes(x = x, y = y, label = add), data = get_sig_data(abund_data = res, diff_data = tmp5), inherit.aes = FALSE) +
	scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
	ylab("Competition index") +
	xlab("") +
	theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 12)) +
	theme(axis.title.y = element_text(size = 18)) +
	theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm")) + 
	labs(color = "Label")

g3



#########################
# Complementarity index
x1 <- genome_interactions$complementarity

res <- data.frame()
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "+", select_matrix = x1, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}
for(i in names(stool_met_network)){
	tmp <- get_matrix_value(network = stool_met_network[[i]], label = "-", select_matrix = x1, group_name = i, module_number = NULL)
	res <- rbind(res, tmp)
}

res$Group %<>% factor(., levels = names(stool_met_network))
res %<>% .[!is.na(.$Value), ]
res$Measure <- res$Group
res$Group <- paste0(res$Measure, "_", res$label)
# Significance
tmp1 <- res
tmp1$Measure <- "Test"
tmp2 <- trans_alpha$new(dataset = NULL)
tmp2$data_alpha <- tmp1
tmp2$group <- "Group"
tmp2$cal_diff(method = "anova", measure = "Test")

tmp3 <- tmp2$res_diff
tmp4 <- tmp3 %>% rownames_to_column(var = "Group")
tmp5 <- res[, c(1, 3, 4)] %>% unique
tmp5 <- left_join(tmp4, tmp5, by = c("Group" = "Group"))
colnames(tmp5)[colnames(tmp5) == "Group"] <- "Type"
colnames(tmp5)[colnames(tmp5) == "label"] <- "Group"
colnames(tmp5)[colnames(tmp5) == "Test"] <- "Significance"

colnames(res)[colnames(res) == "Group"] <- "Type"
colnames(res)[colnames(res) == "label"] <- "Group"
res$Group %<>% factor(., levels = c("+", "-"))
res$Measure %<>% factor(., levels = names(stool_met_network))


g4 <- ggplot(res, aes(x = Measure, y = Value, color = Group)) +
	theme_bw() +
	geom_boxplot(position = position_dodge(0.9)) +
	geom_text(aes(x = x, y = y, label = add), data = get_sig_data(abund_data = res, diff_data = tmp5), inherit.aes = FALSE) +
	scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
	ylab("Complementarity index") +
	xlab("") +
	theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 12)) +
	theme(axis.title.y = element_text(size = 18)) +
	theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm")) + 
	labs(color = "Label")

g4


library(gridExtra)

p1 <- ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2, labels = letters[1:4], font.label = list(size = 20))
cowplot::save_plot("Figure4.png", p1, base_aspect_ratio = 1.6, dpi = 300, base_height = 10)









