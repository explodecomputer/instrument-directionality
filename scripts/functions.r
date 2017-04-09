library(TwoSampleMR)
library(dplyr)


pairwise_directionality <- function(d1, d2, steiger_thresh=0.05)
{
	d1$pval[d1$pval==0] <- 1e-200
	d2$pval[d2$pval==0] <- 1e-200

	l0 <- list()
	for(i in 1:nrow(d1))
	{
		l0[[i]] <- mr_steiger(
			d1$pval[i], 
			d2$pval[i], 
			d1$n[i], 
			d2$n[i]
		)
	}

	dat12 <- recode_dat(data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=d1$bhat,
		beta.outcome=d2$bhat,
		se.exposure=d1$se,
		se.outcome=d2$se,
		pval.exposure=d1$pval,
		pval.outcome=d2$pval,
		samplesize.exposure=d1$n,
		samplesize.outcome=d2$n,
		mr_keep=TRUE
	))

	dat21 <- recode_dat(data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=d2$bhat,
		beta.outcome=d1$bhat,
		se.exposure=d2$se,
		se.outcome=d1$se,
		pval.exposure=d2$pval,
		pval.outcome=d1$pval,
		samplesize.exposure=d2$n,
		samplesize.outcome=d1$n,
		mr_keep=TRUE
	))


	d <- data.frame(
		inst = d1$inst,
		sig1 = d1$pval < 5e-8,
		sig2 = d2$pval < 5e-8,
		dir12 = sapply(l0, function(x) x$correct_causal_direction),
		dir12_p = sapply(l0, function(x) x$steiger_test)
	)

	i12 <- d$sig1 & d$dir12 & d$dir12_p < steiger_thresh
	i21 <- d$sig2 & !d$dir12 & d$dir12_p < steiger_thresh

	res12 <- mr(dat12[i12,], meth="mr_ivw")$pval
	res21 <- mr(dat21[i21,], meth="mr_ivw")$pval

	return(list(d=d, res=c(res12, res21)))
}


find_invalid_instruments <- function(d1, d2, du, steiger_thresh=0.05)
{
	d1$pval[d1$pval==0] <- 1e-200
	d2$pval[d2$pval==0] <- 1e-200
	du$pval[du$pval==0] <- 1e-200

	index <- d1$pval < 5e-8

	# 1. test directionality of each instrument between x and y
	# 2. test directionality of each instrument between x and u
	# 3. using directionality test causality of u on x and x on u
	# 4. if u on x then test directionality of each instrument between u and y
	# 5. test causality of u on y and y on u
	# 6. if u on y and u on x then remove confounding instruments

	xy <- pairwise_directionality(d1, d2, steiger_thresh)
	ux <- pairwise_directionality(du, d1, steiger_thresh)
	uy <- pairwise_directionality(du, d2, steiger_thresh)

	d <- data.frame(
		inst = d1$inst,
		sig_inst = d1$pval < 5e-8,
		xy = xy$d$dir12,
		xy_p = xy$d$dir12_p,
		ux = ux$d$dir12,
		ux_p = ux$d$dir12_p,
		uy = uy$d$dir12,
		uy_p = uy$d$dir12_p,
		ux_sig = ux$res[1] < steiger_thresh,
		uy_sig = uy$res[1] < steiger_thresh
	)

	d$remove_reverse <- 
		d$xy & 
		d$xy_p < steiger_thresh

	d$remove_confounder <- 
		d$ux & 
		d$uy &
		d$ux_p < steiger_thresh & 
		d$uy_p < steiger_thresh &
		ux$res[1] < steiger_thresh &
		uy$res[1] < steiger_thresh
	d$keep <- d$sig_inst & !d$remove_reverse & !d$remove_confounder

	return(d)
}


get_summary_stats <- function(pop1, pop2, popu)
{
	x <- gwas(pop1$x, cbind(pop1$G1, pop1$G2, pop1$Gu))
	y <- gwas(pop2$y, cbind(pop2$G1, pop2$G2, pop2$Gu))
	u <- gwas(popu$u, cbind(popu$G1, popu$G2, popu$Gu))
	x$inst <- y$inst <- u$inst <- rep(c("x", "y", "u"), c(ncol(pop1$G1), ncol(pop1$G2), ncol(pop1$Gu)))

	dat_xy <- recode_dat(data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=x$bhat,
		beta.outcome=y$bhat,
		se.exposure=x$se,
		se.outcome=y$se,
		pval.exposure=x$pval,
		pval.outcome=y$pval,
		samplesize.exposure=x$n,
		samplesize.outcome=y$n,
		mr_keep=TRUE,
		inst = x$inst
	))

	dat_yx <- recode_dat(data.frame(
		exposure="Y",
		id.exposure="Y",
		outcome="X",
		id.outcome="X",
		beta.exposure=y$bhat,
		beta.outcome=x$bhat,
		se.exposure=y$se,
		se.outcome=x$se,
		pval.exposure=y$pval,
		pval.outcome=x$pval,
		samplesize.exposure=y$n,
		samplesize.outcome=x$n,
		mr_keep=TRUE,
		inst = x$inst
	))

	gw <- list(x=x, y=y, u=u)

	return(list(dat_xy=dat_xy, dat_yx=dat_yx, gw=gw))
}


run_sim <- function(nid1, nid2, nidu, ninst1, ninst2, ninstu, var_xy, var_ux, var_uy, var_g1x, var_g2y, var_guu, var_g1y, var_g2x, mu_g1y, mu_g2x, mr_method)
{

	param <- data.frame(nid1 = nid1, nid2 = nid2, nidu = nidu, ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x)

	effs <- make_effs(ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x)

	message("Simulating populations")

	pop1 <- make_pop(effs, nid1)
	pop2 <- make_pop(effs, nid2)
	popu <- make_pop(effs, nidu)

	message("Getting summary stats")

	ss <- get_summary_stats(pop1, pop2, popu)

	message("Selecting instruments")

	xy <- find_invalid_instruments(ss$gw$x, ss$gw$y, ss$gw$u)
	yx <- find_invalid_instruments(ss$gw$y, ss$gw$x, ss$gw$u)

	message("Analysing")

	res <- list()
	for(meth in mr_method)
	{
		message(meth)
		xy_all_res <- with(ss$dat_xy[xy$sig_inst,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
		xy_sel_res <- with(ss$dat_xy[xy$keep,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
		xy_true_res <- with(ss$dat_xy[xy$inst == "x" & xy$sig_inst,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
		yx_all_res <- with(ss$dat_yx[yx$sig_inst,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
		yx_sel_res <- with(ss$dat_yx[yx$keep,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))
		yx_true_res <- with(ss$dat_yx[yx$inst == "y" & yx$sig_inst,], get(meth)(beta.exposure, beta.outcome, se.exposure, se.outcome, default_parameters()))

		res[[meth]] <- data.frame(
			method = meth,
			dir = c("xy", "xy", "xy", "yx", "yx", "yx"),
			inst = c("all", "sel", "true", "all", "sel", "true"),
			b = c(xy_all_res$b, xy_sel_res$b, xy_true_res$b, yx_all_res$b, yx_sel_res$b, yx_true_res$b),
			se = c(xy_all_res$se, xy_sel_res$se, xy_true_res$se, yx_all_res$se, yx_sel_res$se, yx_true_res$se),
			pval = c(xy_all_res$pval, xy_sel_res$pval, xy_true_res$pval, yx_all_res$pval, yx_sel_res$pval, yx_true_res$pval),
			Q = c(xy_all_res$Q, xy_sel_res$Q, xy_true_res$Q, yx_all_res$Q, yx_sel_res$Q, yx_true_res$Q),
			Q_df = c(xy_all_res$Q_df, xy_sel_res$Q_df, xy_true_res$Q_df, yx_all_res$Q_df, yx_sel_res$Q_df, yx_true_res$Q_df),
			nsnp = c(xy_all_res$nsnp, xy_sel_res$nsnp, xy_true_res$nsnp, yx_all_res$nsnp, yx_sel_res$nsnp, yx_true_res$nsnp),
			b_i = NA,
			se_i = NA,
			pval_i = NA
		)

		if(meth == "mr_egger_regression")
		{
			res[[meth]]$b_i <- c(xy_all_res$b_i, xy_sel_res$b_i, xy_true_res$b_i, yx_all_res$b_i, yx_sel_res$b_i, yx_true_res$b_i)
			res[[meth]]$se_i <- c(xy_all_res$se_i, xy_sel_res$se_i, xy_true_res$se_i, yx_all_res$se_i, yx_sel_res$se_i, yx_true_res$se_i)
			res[[meth]]$pval_i <- c(xy_all_res$pval_i, xy_sel_res$pval_i, xy_true_res$pval_i, yx_all_res$pval_i, yx_sel_res$pval_i, yx_true_res$pval_i)
		}
		res[[meth]]$isq <- (res[[meth]]$Q - res[[meth]]$Q_df) / res[[meth]]$Q
	}
	res <- bind_rows(res)

	txy <- table(xy$inst, xy$keep) %>% as.data.frame.matrix()
	tyx <- table(yx$inst, yx$keep) %>% as.data.frame.matrix()
	names(txy) <- names(tyx) <- c("excluded", "included")
	txy$var <- tyx$var <- rownames(txy)
	txy$dir <- "xy"
	tyx$dir <- "yx"
	selection <- rbind(txy, tyx)
	selection$conf_xy <- xy$xu_sig[1] & xy$yu_sig[1]
	selection$conf_yx <- yx$xu_sig[1] & yx$yu_sig[1]

	return(list(selection=selection, res=res, param=param, ss=ss, xy=xy, yx=yx))
}

