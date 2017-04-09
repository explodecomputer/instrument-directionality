library(TwoSampleMR)
library(dplyr)

find_invalid_instruments <- function(d1, d2, d3, steiger_thresh=0.05)
{
	d1$pval[d1$pval==0] <- 1e-200
	d2$pval[d2$pval==0] <- 1e-200
	d3$pval[d3$pval==0] <- 1e-200

	index <- d1$pval < 5e-8

	l0 <- list()
	for(i in 1:nrow(d1))
	{
		l0[[i]] <- mr_steiger(
			d2$pval[i], 
			d1$pval[i], 
			d2$n[i], 
			d1$n[i]
		)
	}
	l1 <- list()
	for(i in 1:nrow(d1))
	{
		l1[[i]] <- mr_steiger(
			d3$pval[i],
			d1$pval[i],
			d3$n[i],
			d1$n[i]
		)
	}
	l2 <- list()
	for(i in 1:nrow(d1))
	{
		l2[[i]] <- mr_steiger(
			d3$pval[i],
			d2$pval[i],
			d3$n[i],
			d2$n[i]
		)
	}

	d <- data.frame(
		inst = d1$inst,
		sig_inst = d1$pval < 5e-8,
		xy = sapply(l0, function(x) x$correct_causal_direction),
		xy_p = sapply(l0, function(x) x$steiger_test),
		ux = sapply(l1, function(x) x$correct_causal_direction),
		ux_p = sapply(l1, function(x) x$steiger_test),
		uy = sapply(l2, function(x) x$correct_causal_direction),
		uy_p = sapply(l2, function(x) x$steiger_test)
	)

	d$remove_reverse <- with(d, xy & xy_p < steiger_thresh)
	d$remove_confounder <- with(d, ux & uy * ux_p < steiger_thresh & uy_p < steiger_thresh)
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

	return(list(selection=selection, res=res, param=param))
}

