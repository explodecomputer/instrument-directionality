library(dplyr)

nom <- paste0("../results/simulate2_", 1:100, ".rdata")
res <- list()
ruck <- list()
validity <- list()
parameters <- list()
param <- list()
for(i in 1:length(nom))
{
	if(file.exists(nom[i]))
	{
		message(i)
		load(nom[i])
		res[[i]] <- lapply(l, function(x) x$res) %>% bind_rows()
		ruck[[i]] <- lapply(l, function(x) x$ruck) %>% bind_rows()
		validity[[i]] <- lapply(l, function(x) x$validity) %>% bind_rows()
		param[[i]] <- lapply(l, function(x) {
			as.data.frame(x$parameters[c("sim", "nidx", "nidy", "nidu", "eff_x.y", "nsnp_x", "nsnp_y", "var_gx.x", "var_gy.y")])
		}) %>% bind_rows()
		parameters[[i]] <- lapply(l, function(x) x$parameters)
	}
}

res <- bind_rows(res)
ruck <- bind_rows(ruck)
validity <- bind_rows(validity)
param <- bind_rows(param)

save(res, ruck, validity, param, parameters, file="../results/simulate2.rdata")
