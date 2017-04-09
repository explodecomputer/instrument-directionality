source("functions.r")
source("simulation-functions.r")

# param <- expand.grid(
# 	ninst1 = 50,
# 	ninst2 = 50,
# 	ninstu = 50,
# 	nid1 = 80000,
# 	nid2 = 80000,
# 	nidu = 80000,
# 	var_xy = c(0.01, 0.1, 0.2),
# 	var_ux = c(0, 0.01, 0.1, 0.2),
# 	var_uy = c(0, 0.01, 0.1, 0.2),
# 	var_g1x = c(0.05, 0.15, 0.3),
# 	var_g2y = c(0.05, 0.15, 0.3),
# 	var_guu = c(0.05, 0.15, 0.3),
# 	var_g1y = c(0, 0.01, 0.05, 0.1),
# 	var_g2x = c(0, 0.01, 0.05, 0.1),
# 	mu_g1y = c(-0.05, 0, 0.05),
# 	mu_g2x = c(-0.05, 0, 0.05)
# )

arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

# chunksize <- ceiling(nrow(param) / chunks)
# t1 <- (jid - 1) * chunksize + 1
# t2 <- min(jid * chunksize, nrow(param))

# message("total size: ", nrow(param))
# message("running: ", t1, " to ", t2)

# param <- param[t1:t2, ]


message("running ", chunks, " simulations")
set.seed(jid)

l <- list()
for(i in 1:chunks)
{
	message(i)
	# Sizes
	ninst1 <- 50
	ninst2 <- 50
	ninstu <- 50
	nid1 <- 80000
	nid2 <- 80000
	nidu <- 80000

	# Causal effects
	var_xy <- runif(1, 0.01, 0.2)
	var_ux <- runif(1, 0.01, 0.2)
	var_uy <- runif(1, 0.01, 0.2)

	# Genetic effects
	var_g1x <- runif(1, 0.05, 0.3)
	var_g2y <- runif(1, 0.05, 0.3)
	var_guu <- runif(1, 0.05, 0.3)

	# Horizontal pleiotropy
	# var_g1y <- min(rbeta(1, 1, 100), 0.2)
	var_g1y <- 0
	# var_g2x <- min(rbeta(1, 1, 100), 0.2)
	var_g2x <- 0
	# mu_g1y <- min(rbeta(1, 1, 200), 0.05) * sample(c(1,-1), 1)
	# mu_g2x <- min(rbeta(1, 1, 300), 0.05) * sample(c(1,-1), 1)
	mu_g1y <- 0
	mu_g2x <- 0

	mr_method <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")

	l[[i]] <- run_sim(nid1 = nid1, nid2 = nid2, nidu = nidu, ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x, mr_method = mr_method)
	l[[i]]$selection$sim <- i + chunks * (jid - 1)
	l[[i]]$res$sim <- i + chunks * (jid - 1)
	l[[i]]$param$sim <- i + chunks * (jid - 1)
}

save(l, file=out)