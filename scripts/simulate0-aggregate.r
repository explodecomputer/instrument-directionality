library(dplyr)

nom <- paste0("../results/simulate0_", 1:100, ".rdata")
m1 <- list()
m2 <- list()
m3 <- list()
m4 <- list()
m5 <- list()
for(i in 1:length(nom))
{
	message(i)
	load(nom[i])

	param <- bind_rows(lapply(l, function(x) x$param))
	param$rucker <- sapply(l, function(x) x$res)
	rucker <- bind_rows(lapply(l, function(x) x$rucker))
	intercept <- bind_rows(lapply(l, function(x) x$intercept))
	Q <- bind_rows(lapply(l, function(x) x$Q))
	mrmode <- bind_rows(lapply(l, function(x) x$mode))
	mrmedian <- bind_rows(lapply(l, function(x) x$median))

	index <- param$sim
	rucker$sim <- rep(index, each=4)
	intercept$sim <- rep(index, each=2)
	Q$sim <- rep(index, each=3)
	mrmode$sim <- rep(index, each=4)
	mrmedian$sim <- rep(index, each=3)

	m1[[i]] <- merge(param, rucker, by="sim")
	m2[[i]] <- merge(param, intercept, by="sim")
	m3[[i]] <- merge(param, Q, by="sim")
	m4[[i]] <- merge(param, mrmedian, by="sim")
	m5[[i]] <- merge(param, mrmode, by="sim")
}

rucker <- bind_rows(m1)
rucker_sel <- group_by(rucker, sim) %>%
	do({
		.[LETTERS[1:4] %in% .$rucker, ]
	})
intercept <- bind_rows(m2)
Q <- bind_rows(m3)
mrmedian <- bind_rows(m4)
mrmode <- bind_rows(m5)

res <- rbind(
	data.frame
)

save(dat, dat_sel, intercept, Q, mrmedian, mrmode, file="../results/simulate0.rdata")
