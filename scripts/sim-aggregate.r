library(dplyr)

nom <- paste0("../results/simulate3_", 1:500, "-metrics.rdata")
res <- list()
metrics <- list()
validity <- list()
parameters <- list()
param <- list()
for(i in 1:length(nom))
{
	if(file.exists(nom[i]))
	{
		message(i)
		load(nom[i])
		res[[i]] <- lapply(m, function(x) x$res) %>% bind_rows()
		metrics[[i]] <- lapply(m, function(x) x$metrics) %>% bind_rows()
		validity[[i]] <- lapply(m, function(x) x$validity) %>% bind_rows()
		param[[i]] <- lapply(m, function(x) {
			as.data.frame(x$parameters[c("sim", "nidx", "nidy", "eff_x.y", "nsnp_x", "nsnp_y", "var_gx.x", "var_gy.y", "var_gx.y", "var_gy.x", "prop_gx.y", "prop_gy.x", "mu_gx.y", "mu_gy.x")])
		}) %>% bind_rows()
		parameters[[i]] <- lapply(m, function(x) x$parameters)
	}
}

res <- bind_rows(res)
metrics <- bind_rows(metrics)
validity <- bind_rows(validity)
param <- bind_rows(param)

save(res, metrics, validity, param, parameters, file="../results/simulate3.rdata")

