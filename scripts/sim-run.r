library(TwoSampleMR)
library(dplyr)
source("fun-simulations.r")


arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
sims <- as.numeric(arguments[2])
out <- arguments[3]

message("running ", sims, " simulations")
set.seed(jid)

l <- list()
j <- 1
for(i in 1:sims)
{
	k <- (jid - 1) * sims + i
	message(i, " (", k, ")")
	ss <- try(create_system(
		nidx=sample(20000:500000, 1), 
		nidy=sample(20000:500000, 1), 
		nidu=0, 
		nu=sample(0:10, 1), 
		na=0,
		nb=0,
		var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
		nsnp_x=sample(1:200, 1),
		nsnp_y=sample(1:200, 1),
		var_gx.x=runif(1, 0.01, 0.1),
		var_gy.y=runif(1, 0.01, 0.1),
		var_gx.y=runif(1, 0.001, 0.01),
		mu_gx.y=runif(1, -0.005, 0.005),
		prop_gx.y=runif(1, 0, 1),
		var_gy.x=runif(1, 0.001, 0.01),
		mu_gy.x=runif(1, -0.005, 0.005),
		prop_gy.x=runif(1, 0, 1)
	))
	if(class(ss) != "try-error")
	{
		l[[j]] <- ss
		j <- j + 1
	}
}

save(l, file=out)
