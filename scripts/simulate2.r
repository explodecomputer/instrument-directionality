library(TwoSampleMR)
library(dplyr)
source("simulation-functions.r")
source("strategies.r")


arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

message("running ", chunks, " simulations")
set.seed(jid)



test_ss <- function(ss, i)
{
	ss1 <- strategies_oracle(ss)
	ss2 <- strategies_tophits(ss)
	ss3 <- strategies_directionality(ss)

	res1xy <- mr_all(ss1$xy)
	res1yx <- mr_all(ss1$yx)
	res2xy <- mr_all(ss2$xy)
	res2yx <- mr_all(ss2$yx)
	res3xy <- mr_all(ss3$xy)
	res3yx <- mr_all(ss3$yx)

	# Estimates
	res <- bind_cols(
		bind_rows(res1xy$res, res1yx$res, res2xy$res, res2yx$res, res3xy$res, res3yx$res),
		expand.grid(1:nrow(res1xy$res), hypothesis=c("xy", "yx"), strategy=c("oracle", "tophits", "steiger"))[,-1]
	)
	res$sim <- i

	# Rucker models
	ruck <- bind_cols(
		data.frame(model = c(res1xy$rucker$res, res1xy$ruckercd$res, res1yx$rucker$res, res1yx$ruckercd$res, res2xy$rucker$res, res2xy$ruckercd$res, res2yx$rucker$res, res2yx$ruckercd$res, res3xy$rucker$res, res3xy$ruckercd$res, res3yx$rucker$res, res3yx$ruckercd$res)),
		expand.grid(method=c("Rucker", "Rucker (No outliers)"), hypothesis=c("xy", "yx"), strategy=c("oracle", "tophits", "steiger"))
	)
	ruck$sim <- i

	# Parameters
	p <- ss$parameters
	p$nidx <- ss$x$x$n[1]
	p$nidy <- ss$y$y$n[1]
	p$nidu <- ss$u[[2]]$x$n[1]
	p$sim <- i

	return(list(res=res, ruck=ruck, parameters=p))
}


l <- list()
for(i in 1:chunks)
{
	j <- (jid - 1) * chunks + i
	message(i, " (", j, ")")
	ss <- try(create_system(
		nidx=sample(20000:100000, 1), 
		nidy=sample(20000:100000, 1), 
		nidu=sample(2000:10000, 1), 
		nu=5, na=5, nb=5, 
		var_gx.x=runif(1, 0.01, 0.1),
		var_gy.y=runif(1, 0.01, 0.1),
		var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1)
	))
	l[[i]] <- try(test_ss(ss, j))
}

save(l, file=out)
