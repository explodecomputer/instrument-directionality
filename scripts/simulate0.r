library(TwoSampleMR)
source("simulation-functions.r")

param <- expand.grid(
	nsim = c(1:100),
	nsnp = c(5,20,50),
	nid1 = 50000,
	nid2 = 50000,
	var_xy = c(0, 0, 0, 0, 0, 0.005, 0.01, 0.025, 0.05, 0.1),
	var_g1x = c(0.001, 0.05, 0.1),
	var_g1y = c(0, 0.001, 0.025),
	mu_g1y = c(-0.05, 0, 0.05)
)
param$sim <- 1:nrow(param)

arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]
set.seed(jid)

chunksize <- ceiling(nrow(param) / chunks)
t1 <- (jid - 1) * chunksize + 1
t2 <- min(jid * chunksize, nrow(param))

message("total size: ", nrow(param))
message("running: ", t1, " to ", t2)

param <- param[t1:t2, ]


l <- list()
for(i in 1:nrow(param))
{
	message(i)
	effs <- make_effs(ninst1=param$nsnp[i], var_xy=param$var_xy[i], var_g1x=param$var_g1x[i], mu_g1y=param$mu_g1y[i])
	pop1 <- make_pop(effs, param$nid1[i])
	pop2 <- make_pop(effs, param$nid2[i])
	dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
	dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
	dat <- recode_dat(make_dat(dat1, dat2))
	res <- mr_all(dat)
	res$res$sim <- param$sim[i]
	res$param <- param[i,]
	res$dat <- dat
	l[[i]] <- res
}

save(l, file=out)
