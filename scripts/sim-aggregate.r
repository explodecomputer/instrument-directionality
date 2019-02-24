library(dplyr)

args <- commandArgs(T)
input_dir <- args[1]

input <- file.path(input_dir, list.files(input_dir, pattern="sim_*")) %>% grep("sim_", ., value=TRUE)

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


pleiotropy_metric <- function(x)
{
	a <- x$parameters
	if(is.null(a)) return(NULL)

	a$id <- x$id
	instx <- data_frame(
		id = paste0(a$id, "x"),
		SNP = paste0("i", 1:length(a$eff_gx.x)),
		eff_d = a$eff_gx.x^2,
		eff_p = a$eff_gx.y^2
	)
	revx <- data_frame(
		id = paste0(a$id, "x"),
		SNP = paste0("r", 1:length(a$eff_gy.x)),
		eff_d = a$eff_gy.x^2,
		eff_p = a$eff_gy.y^2
	)
	insty <- data_frame(
		id = paste0(a$id, "y"),
		SNP = paste0("i", 1:length(a$eff_gy.y)),
		eff_d = a$eff_gy.y^2,
		eff_p = a$eff_gy.x^2
	)
	revy <- data_frame(
		id = paste0(a$id, "y"),
		SNP = paste0("r", 1:length(a$eff_gx.y)),
		eff_d = a$eff_gx.y^2 + a$eff_gx.x^2 * a$eff_x.y^2,
		eff_p = a$eff_gx.x^2
	)
	o <- bind_rows(instx, insty, revx, revy)
	if(length(a$u) > 0)
	{
		confx <- data_frame(
			id = paste0(a$id, "x"),
			eff_d = lapply(a$u, function(x) (x$eff_u.x * x$eff_gu.u)^2) %>% unlist,
			eff_p = lapply(a$u, function(x) (x$eff_u.y * x$eff_gu.u)^2) %>% unlist,
		)
		confx$SNP <- paste0("u", 1:nrow(confx))
		confy <- data_frame(
			id = paste0(a$id, "y"),
			eff_d = lapply(a$u, function(x) (x$eff_u.y * x$eff_gu.u)^2) %>% unlist,
			eff_p = lapply(a$u, function(x) (x$eff_u.x * x$eff_gu.u)^2) %>% unlist,
		)
		confy$SNP <- paste0("u", 1:nrow(confy))
		o <- bind_rows(o, confx, confy)
	}
	o <- group_by(o, id, w=substr(SNP, 1,1)) %>%
	summarise(
		nonpl = sum((eff_d / (eff_d + eff_p)) * (eff_d / sum(eff_d))),
		effd = sum(eff_d)
	)

	return(o)
}


estimates <- list()
info <- list()
instrument_validity <- list()
heterogeneity <- list()
directional_pleiotropy <- list()
param <- list()
parameters <- list()
plei_summary <- list()

# mclapply(1:length(input), function(i)
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

		parameters[sapply(parameters, is.null)] <- NULL
		parameters <- lapply(parameters, function(x) {
			x[sapply(x, is.null)] <- NULL
			return(x)
		})

		plei_summary[[i]] <- lapply(l, function(x) pleiotropy_metric(x)) %>% bind_rows
		parameters[[i]] <- lapply(l, function(x) {
			x$parameters$id <- x$id
			x$parameters
		})
	}
}

estimates <- bind_rows(estimates)
info <- bind_rows(info)
instrument_validity <- bind_rows(instrument_validity)
heterogeneity <- bind_rows(heterogeneity)
directional_pleiotropy <- bind_rows(directional_pleiotropy)
param <- bind_rows(param)
plei_summary <- bind_rows(plei_summary)

save(estimates, file=file.path(input_dir, "agg-estimates.rdata"))
save(info, file=file.path(input_dir, "agg-info.rdata"))
save(instrument_validity, file=file.path(input_dir, "agg-instrument_validity.rdata"))
save(heterogeneity, file=file.path(input_dir, "agg-heterogeneity.rdata"))
save(directional_pleiotropy, file=file.path(input_dir, "agg-directional_pleiotropy.rdata"))
save(param, file=file.path(input_dir, "agg-param.rdata"))
save(plei_summary, file=file.path(input_dir, "agg-plei_summary.rdata"))



subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered) %>% summarise(ncorrect=sum(beta_correct, na.rm=T), n=n())

subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered) %>% summarise(ncorrect=sum(beta_best, na.rm=T), n=n())

subset(estimates, selection != "oracle") %>% group_by(steiger_filtered, outlier_filtered, method) %>% summarise(ncorrect=sum(beta_correct, na.rm=T), n=n(), coef=mean(as.numeric(beta_correct) / se,na.rm=T)) %>% arrange(desc(coef)) %>% as.data.frame


