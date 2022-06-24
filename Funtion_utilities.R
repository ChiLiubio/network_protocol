
library(magrittr)
# network: must be trans_network object
# select_matrix: matrix used for the value extraction; colnames should be same with rownames
# group_name: the identifier for the result
# module_number: how many modules are used; default NULL
get_matrix_value <- function(network, label = c("+", "-"), select_matrix, group_name, module_number = NULL){
	if(is.null(module_number)){
		tmp <- network$res_edge_table %>% .[.$label %in% label, ]
		if(nrow(tmp) == 0){
			res <- NA
		}else{
			select_value <- lapply(seq_len(nrow(tmp)), function(x){
				if(all(c(tmp[x, 1], tmp[x, 2]) %in% colnames(select_matrix))){
					select_matrix[tmp[x, 1], tmp[x, 2]]
				}else{
					NA
				}
			}) %>% unlist
			res <- data.frame(Group = group_name, Value = select_value, label = paste0(label, collapse = "|"))
		}
	}else{
		if(!is.numeric(module_number)){
			stop("module_number must be numeric!")
		}
		if(is.null(network$res_node_table)){
			network$get_node_table(node_roles = FALSE)
		}
		if(! "module" %in% colnames(network$res_node_table)){
			stop("please first use cal_module function to calculate modularity!")
		}
		# check min module number
		if(length(unique(network$res_node_table$module)) < module_number){
			module_number <- length(unique(network$res_node_table$module))
		}
		use_modules <- paste0("M", 1:module_number)
		res <- NULL
		for(k in use_modules){
			module_nodes <- network$res_node_table %>% .[.$module == k, ] %>% rownames
			t1 <- clone(network)
			t1$res_network <- t1$subset_network(node = module_nodes, rm_single = TRUE)
			suppressMessages(t1$get_edge_table())
			tmp <- t1$res_edge_table %>% .[.$label %in% label, ]
			if(nrow(tmp) == 0){
				next
			}else{
				select_value <- lapply(seq_len(nrow(tmp)), function(x){
					if(all(c(tmp[x, 1], tmp[x, 2]) %in% colnames(select_matrix))){
						select_matrix[tmp[x, 1], tmp[x, 2]]
					}else{
						NA
					}
				}) %>% unlist
				res <- rbind(res, data.frame(Group = group_name, Value = select_value, module = k, label = paste0(label, collapse = "|")))
			}
		}
	}
	res
}


get_sig_data <- function(abund_data, diff_data){
	add_sig_label = "Significance"
	barwidth = 0.9
	x_axis_order <- levels(abund_data$Group)

	x_mid <- c()
	annotations <- c()
	y_position <- c()

	start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
	increase_bar_mid <- barwidth/length(x_axis_order)
	all_taxa <- levels(abund_data$Measure)

	for(j in all_taxa){
		select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Measure == j, ]
		for(i in seq_len(nrow(select_use_diff_data))){
			# first determine the bar range
			mid_num <- match(j, all_taxa) - 1
			annotations %<>% c(., select_use_diff_data[i, add_sig_label])
			# middle if only one group found
			if(nrow(select_use_diff_data) == 1){
				x_mid %<>% c(., mid_num + start_bar_mid + 0.5 * increase_bar_mid)
			}else{
				x_mid %<>% c(., mid_num + (start_bar_mid + (match(select_use_diff_data[i, "Group"], x_axis_order) - 1) * increase_bar_mid))
			}
			abund_data_select <- abund_data[abund_data$Group == select_use_diff_data[i, "Group"] & abund_data$Measure == j, ]
			y_position %<>% c(., 1.05 * max(abund_data_select$Value))						
		}
	}
	textdf <- data.frame(
		x = x_mid, 
		y = y_position, 
		add = annotations, 
		stringsAsFactors = FALSE
		)
	textdf
}






