## -- coefplot fun -- ## 

coef_plot = function(mod, x_se = 2) {
	
	variable_names = rownames(mod$stderr.coefficients)
	mod$stderr.coefficients %>% 
		as.data.frame() %>% 
		mutate(variable = variable_names) %>% 
		pivot_longer(cols = !(variable)) -> 
		SE
	mod$coefficients%>% 
		as.data.frame() %>% 
		mutate(variable = variable_names) %>% 
		pivot_longer(cols = !(variable), values_to = "estimate", 
		             names_to = "taxa") -> 
		CO
	CO$SE = SE$value
	CO %<>% 
		mutate(lb = estimate - x_se*SE, 
		       ub = estimate + x_se*SE)
	CO %<>% 
		mutate(significant_lb = sign(estimate) == sign(lb), 
		       significant_ub = sign(estimate) == sign(ub),
		       significant = significant_lb & significant_ub)
	
	# If not significant set estimate to zero 
	CO %<>% mutate(
		estimate = ifelse(significant, estimate, 0)
	)
	
	# convert back ot wide format 
	CO %>% 
		dplyr::select(!c("SE", "lb","ub","significant_lb" ,"significant_ub", "significant")) %>% 
		pivot_wider(id_cols = variable, 
		            names_from = taxa,
		            values_from = estimate) %>% 
		as.data.frame() -> 
		coef_plot 
	rownames(coef_plot) = coef_plot$variable
	coef_plot %<>% filter(variable != "(Intercept)")
	coef_plot %<>% dplyr::select(-variable)
	
	
	while (any (colSums(coef_plot) == 0)) {
		col_id = which(colSums(coef_plot) == 0)
		if (length(col_id)>0)
			coef_plot = coef_plot[, -col_id]
		row_id = which(rowSums(coef_plot) == 0)
		if (length(row_id)>0)
			coef_plot = coef_plot[-row_id, ]
	}
	return(coef_plot)
}