library(tidyverse)
library(TwoSampleMR)

Isq <- function(y,s)
{
	k <- length(y)
	w <- 1/s^2
	sum.w <- sum(w)
	mu.hat <- sum(y*w)/sum.w
	Q <- sum(w*(y-mu.hat)^2)
	Isq <- (Q - (k-1))/Q
	Isq <- max(0,Isq)
	return(Isq)
}

instrument_validity <- function(ssn)
{
	data.frame(
		hypothesis=rep(c("xy", "yx"), each=3),
		type=rep(c("valid", "confounder", "reverse"), times=2),
		value=c(
			sum(ssn$xy$inst == "x"),
			sum(grepl("u", ssn$xy$inst)),
			sum(ssn$xy$inst == "y"),
			sum(ssn$yx$inst == "y"),
			sum(grepl("u", ssn$yx$inst)),
			sum(ssn$yx$inst == "x")
		)
	)
}

best_model <- function(resall, bxy)
{
	res <- subset(resall, strategy != "oracle")
	res_xy <- res[res$hypothesis == "xy", ]
	res_yx <- res[res$hypothesis == "yx", ]

	res_xy$beta_correct <- res_xy$CI_low <= bxy & res_xy$CI_upp >= bxy
	res_xy$beta_best <- FALSE
	res_xy$beta_best[which.min(abs(res_xy$Estimate - bxy))] <- TRUE
	res_xy$pval_sig <- res_xy$P < 0.05
	res_xy$pval_lowest <- FALSE
	res_xy$pval_lowest[which.min(res_xy$P)] <- TRUE
	res_xy$pval_highest <- FALSE
	res_xy$pval_highest[which.max(res_xy$P)] <- TRUE

	res_yx$beta_correct <- res_yx$CI_low <= 0 & res_yx$CI_upp >= 0
	res_yx$beta_best <- FALSE
	res_yx$beta_best[which.min(abs(res_yx$Estimate - 0))] <- TRUE
	res_yx$pval_sig <- res_yx$P < 0.05
	res_yx$pval_lowest <- FALSE
	res_yx$pval_lowest[which.min(res_yx$P)] <- TRUE
	res_yx$pval_highest <- FALSE
	res_yx$pval_highest[which.max(res_yx$P)] <- TRUE
	res <- rbind(res_xy, res_yx, subset(resall, strategy == "oracle"))
	return(res)
}

test_ss <- function(ss, i)
{
	ss1 <- strategies_oracle(ss)
	ss2 <- strategies_tophits(ss)
	ss3 <- strategies_directionality(ss)

	metrics1 <- get_metrics(ss2[[1]]) %>% as.data.frame()
	metrics1$hypothesis <- "xy"
	metrics1$sim <- i
	metrics2 <- get_metrics(ss2[[2]]) %>% as.data.frame()
	metrics2$hypothesis <- "yx"
	metrics2$sim <- i
	metrics <- bind_rows(metrics1, metrics2)

	res1xy <- run_mr(ss1$xy)
	res1yx <- run_mr(ss1$yx)
	res2xy <- run_mr(ss2$xy)
	res2yx <- run_mr(ss2$yx)
	res3xy <- run_mr(ss3$xy)
	res3yx <- run_mr(ss3$yx)

	# Estimates
	res <- bind_cols(
		bind_rows(res1xy, res1yx, res2xy, res2yx, res3xy, res3yx),
		expand.grid(1:nrow(res1xy), hypothesis=c("xy", "yx"), strategy=c("oracle", "tophits", "steiger"))[,-1]
	)
	res$sim <- i

	# Parameters
	p <- ss$parameters
	p$nidx <- ss$x$x$n[1]
	p$nidy <- ss$y$y$n[1]
	p$nidu <- ss$u[[2]]$x$n[1]
	p$sim <- i

	# Instrument validity
	v1 <- instrument_validity(ss1)
	v2 <- instrument_validity(ss2)
	v3 <- instrument_validity(ss3)
	validity <- bind_cols(
		bind_rows(v1, v2, v3),
		expand.grid(1:nrow(v1), strategy=c("oracle", "tophits", "steiger"))
	)
	validity$sim <- i
	validity <- subset(validity, select=-c(Var1))

	res <- best_model(res, p$eff_x.y)
	return(list(res=res, parameters=p, validity=validity, metrics=metrics))
}

source("strategies.r")
source("simulation-functions.r")

###

infile <- commandArgs(T)[1]
outfile <- commandArgs(T)[2]
jid <- as.numeric(commandArgs(T)[3])

load(infile)

m <- list()
sims <- length(l)
for(i in 1:sims)
{
	message(i)
	k <- (jid - 1) * sims + i
	a <- try(test_ss(l[[i]], k))
	if(class(a) != "try-error")
	{
		m[[i]] <- a
	}
}

save(m, l, file=outfile)
