library(dplyr)

nom <- paste0("../results/simulate1_", 1:100, ".rdata")
m1 <- list()
m2 <- list()
for(i in 1:length(nom))
{
	message(i)
	load(nom[i])

	m1[[i]] <- lapply(l, function(x)
	{
		merge(x$param, x$res, by="sim")
	}) %>% bind_rows()

	m2[[i]] <- lapply(l, function(x)
	{
		merge(x$param, x$selection, by="sim")
	}) %>% bind_rows()
}

res <- bind_rows(m1)
selection <- bind_rows(m2)

save(res, selection, file="../results/simulate1.rdata")
