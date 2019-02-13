library(dplyr)

args <- commandArgs(T)
input_dir <- args[1]
output <- args[2]

input <- file.path(input_dir, list.files(input_dir, pattern="*metrics.rdata"))

estimates <- list()
info <- list()
instrument_validity <- list()
heterogeneity <- list()
directional_pleiotropy <- list()
param <- list()

check_nulls <- function(l, node)
{
	l <- lapply(l, function(x){
		if(class(x$estimates) == "try-error") x$estimates <- list(ox=list(),oy=list(),ex=list(),ey=list())
		if(is.null(x$estimates$ox[[node]])) x$estimates$ox[[node]] <- data_frame()
		if(is.null(x$estimates$oy[[node]])) x$estimates$oy[[node]] <- data_frame()
		if(is.null(x$estimates$ex[[node]])) x$estimates$ex[[node]] <- data_frame()
		if(is.null(x$estimates$ey[[node]])) x$estimates$ey[[node]] <- data_frame()
		return(x)
	})
	return(l)
}

aggregate <- function(l, node)
{
	out <- lapply(l, function(x)
		bind_rows(
			x$estimates$ox[[node]] %>% mutate(hypothesis = "x", selection = "o"),
			x$estimates$oy[[node]] %>% mutate(hypothesis = "y", selection = "o"),
			x$estimates$ex[[node]] %>% mutate(hypothesis = "x", selection = "e"),
			x$estimates$ey[[node]] %>% mutate(hypothesis = "y", selection = "e")
		) %>% mutate(id = x$id)
	) %>% bind_rows
	return(out)
}

for(i in 1:length(input))
{
	if(file.exists(input[i]))
	{
		message(i)
		load(input[i])

		l <- check_nulls(l, "directional_pleiotropy")
		l <- check_nulls(l, "heterogeneity")
		l <- check_nulls(l, "estimates")
		l <- check_nulls(l, "info")

		estimates[[i]] <- aggregate(l, "estimates")
		heterogeneity[[i]] <- aggregate(l, "heterogeneity")
		info[[i]] <- aggregate(l, "info")
		directional_pleiotropy[[i]] <- aggregate(l, "directional_pleiotropy")

		instrument_validity[[i]] <- lapply(l, function(x) x$estimates$instrument_validity) %>% bind_rows

		param[[i]] <- lapply(l, function(x) {
			as_data_frame(x$parameters[c("eff_x.y", "nsnp_x", "nsnp_y", "var_gx.x", "var_gy.y", "var_gx.y", "var_gy.x", "prop_gx.y", "prop_gy.x", "mu_gx.y", "mu_gy.x")]) %>% 
			mutate(id = x$id, nidx = x$x$x$n[1], nidy = x$y$y$n[1])
		}) %>% bind_rows()
	}
}

estimates <- bind_rows(estimates)
info <- bind_rows(info)
instrument_validity <- bind_rows(instrument_validity)
heterogeneity <- bind_rows(heterogeneity)
directional_pleiotropy <- bind_rows(directional_pleiotropy)
param <- bind_rows(param)

save(estimates, info, instrument_validity, heterogeneity, directional_pleiotropy, param, file=output)




subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered) %>% summarise(ncorrect=sum(beta_correct, na.rm=T), n=n())

subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered) %>% summarise(ncorrect=sum(beta_best, na.rm=T), n=n())

subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered, method) %>% summarise(ncorrect=sum(beta_correct, na.rm=T), n=n(), coef=mean(as.numeric(beta_correct) / se,na.rm=T)) %>% arrange(desc(coef)) %>% as.data.frame


