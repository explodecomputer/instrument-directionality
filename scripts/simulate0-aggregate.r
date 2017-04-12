library(dplyr)

nom <- paste0("../results/simulate0_", 1:100, ".rdata")
chunksize <- 108
m1 <- list()
m2 <- list()
m3 <- list()
for(i in 1:length(nom))
{
	message(i)
	load(nom[i])

	param <- bind_rows(lapply(l, function(x) x$param))
	param$rucker <- sapply(l, function(x) x$res)
	dat <- bind_rows(lapply(l, function(x) x$dat))
	intercept <- bind_rows(lapply(l, function(x) x$intercept))
	Q <- bind_rows(lapply(l, function(x) x$Q))

	index <- chunksize * (i - 1) + 1:chunksize
	param$sim <- index
	dat$sim <- rep(index, each=4)
	intercept$sim <- rep(index, each=2)
	Q$sim <- rep(index, each=3)

	m1[[i]] <- merge(param, dat, by="sim")
	m2[[i]] <- merge(param, intercept, by="sim")
	m3[[i]] <- merge(param, Q, by="sim")
}

dat <- bind_rows(m1)
dat_sel <- group_by(dat, sim) %>%
	do({
		.[LETTERS[1:4] %in% .$rucker, ]
	})
intercept <- bind_rows(m2)
Q <- bind_rows(m3)

save(dat, dat_sel, intercept, Q, file="../results/simulate0.rdata")
