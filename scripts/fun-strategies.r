
steiger_simple <- function(p_exp, p_out, n_exp, n_out)
{
	p_exp[p_exp == 0] <- 1e-200
	p_out[p_out == 0] <- 1e-200
	r_exp <- sqrt(sum(get_r_from_pn(p_exp, n_exp)^2))
	r_out <- sqrt(sum(get_r_from_pn(p_out, n_out)^2))
	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	return(list(dir=r_exp > r_out, pval=rtest$p))
}

# Strategies
#
# - Oracle - only true strong instruments
# - Top hits
# - Remove x-y directionality
# - Scan through confounders, identify new instruments evaluate directionality of u-x, evaluate if u causes x, if u causes x then do MR of all u instruments' effects on X against Y


# This simulation 
strategies_oracle <- function(ss)
{
	# For x on y
	d1 <- ss$x$x
	d2 <- ss$y$y

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

	dat12$inst <- d1$inst
	dat21$inst <- d2$inst

	p <- sapply(ss$parameters$u, function(x) x$name_u)
	xu <- which(p == "a")
	yu <- which(p == "b")

	dat12 <- subset(dat12, pval.exposure < 5e-8 & (inst == "x" | inst %in% paste0("u", xu)))
	dat21 <- subset(dat21, pval.exposure < 5e-8 & (inst == "y" | inst %in% paste0("u", yu)))

	return(list(xy=dat12, yx=dat21))

}

strategies_tophits <- function(ss)
{
	# For x on y
	d1 <- ss$x$x
	d2 <- ss$y$y

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

	dat12$inst <- d1$inst
	dat21$inst <- d2$inst

	dat12 <- subset(dat12, pval.exposure < 5e-8)
	dat21 <- subset(dat21, pval.exposure < 5e-8)

	return(list(xy=dat12, yx=dat21))
}

strategies_directionality <- function(ss, steiger_thresh=0.05)
{
	# For x on y
	d1 <- ss$x$x
	d2 <- ss$y$y
	d1$pval[d1$pval==0] <- 1e-200
	d2$pval[d2$pval==0] <- 1e-200

	l0 <- list()
	for(i in 1:nrow(d1))
	{
		l0[[i]] <- steiger_simple(
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


	dat12$inst <- d1$inst
	dat12$sig <- d1$pval < 5e-8
	dat12$dir <- sapply(l0, function(x) x$dir)
	dat12$dir_p <- sapply(l0, function(x) x$pval)
	dat12$mr_keep <- dat12$sig & dat12$dir & dat12$dir_p < steiger_thresh
	dat12 <- subset(dat12, mr_keep)

	dat21$inst <- d2$inst
	dat21$sig <- d2$pval < 5e-8
	dat21$dir <- ! sapply(l0, function(x) x$dir)
	dat21$dir_p <- sapply(l0, function(x) x$pval)
	dat21$mr_keep <- dat21$sig & dat21$dir & dat21$dir_p < steiger_thresh
	dat21 <- subset(dat21, mr_keep)

	return(list(xy=dat12, yx=dat21))
}

# Scan through confounders, identify new instruments evaluate directionality of u-x, evaluate if u causes x, if u causes x then do MR of all u instruments' effects on X against Y
# Still needs work
strategies_spider <- function(ss, steiger_thresh=0.05)
{
	# For x on y
	d1 <- ss$x$x
	d2 <- ss$y$y
	d1$pval[d1$pval==0] <- 1e-200
	d2$pval[d2$pval==0] <- 1e-200

	l0 <- list()
	for(i in 1:nrow(d1))
	{
		l0[[i]] <- steiger_simple(
			d1$pval[i], 
			d2$pval[i], 
			d1$n[i], 
			d2$n[i]
		)
	}

	d1 <- subset(d1, sapply(l0, function(x) x$dir) & sapply(l0, function(x) x$pval < steiger_thresh) & pval < 5e-8)
	d2 <- subset(d2, ! sapply(l0, function(x) x$dir) & sapply(l0, function(x) x$pval < steiger_thresh) & pval < 5e-8)

	# For each instrument find if it associates with any other a/b/u
	d1$spider <- NA
	d1$rsq <- get_r_from_pn(d1$pval, d1$n)^2
	m <- list()
	for(i in 1:nrow(d1))
	{
		# Look through Us
		if(length(ss$u) > 0)
		{
			lu <- ss$u[[1]][grepl("u", names(ss$u[[1]]))]
			lu <- lapply(lu, function(x)
			{
				x <- subset(x, inst == d1$inst[i] & snp == d1$snp[i])
				# x$rsq <- 
			})
		}

		# Look through As
		if(length(ss$a) > 0)
		{
			la <- ss$a[[1]][grepl("a", names(ss$a[[1]]))]
			la <- lapply(la, function(x)
			{
				subset(x, inst == d1$inst[i] & snp == d1$snp[i])
			})
		}

		# Look through Bs
		if(length(ss$b) > 0)
		{
			lb <- ss$b[[1]][grepl("b", names(ss$b[[1]]))]
			lb <- lapply(lb, function(x)
			{
				subset(x, inst == d1$inst[i] & snp == d1$snp[i])
			})
		}
		l <- lapply(c(lu, la, lb), function(x) {
			subset(x)
		})
	}
}


system_metrics <- function(dat)
{
	library(car)

	# Number of SNPs
	# Sample size outcome
	# Sample size exposure
	metrics <- list()
	metrics$nsnp <- nrow(dat)
	metrics$nout <- mean(dat$samplesize.outcome, na.rm=TRUE)
	metrics$nexp <- mean(dat$samplesize.exposure, na.rm=TRUE)

	# F stats
	Fstat <- qf(dat$pval.exposure, 1, dat$samplesize.exposure, lower.tail=FALSE)
	Fstat[is.infinite(Fstat)] <- 300
	metrics$meanF <- mean(Fstat, na.rm=TRUE)
	metrics$varF <- var(Fstat, na.rm=TRUE)
	metrics$medianF <- median(Fstat, na.rm=TRUE)

	# IF more than 1 SNP

	if(nrow(dat) > 1)
	{
		# Egger-Isq
		metrics$egger_isq <- Isq(dat$beta.exposure, dat$se.exposure)
	}

	if(nrow(dat) > 2)
	{	
		# IF more than 2 SNP
		ruck <- mr_rucker(dat)

		# Q_ivw
		# Q_egger
		# Q_diff
		metrics$Isq <- (ruck$Q$Q[1] - (ruck$Q$df[1]-1))/ruck$Q$Q[1]
		metrics$Isqe <- (ruck$Q$Q[2] - (ruck$Q$df[2]-1))/ruck$Q$Q[2]
		metrics$Qdiff <- ruck$Q$Q[3]

		# Intercept / se
		metrics$intercept <- abs(ruck$intercept$Estimate[1]) / ruck$intercept$SE[1]

		# Influential outliers
		dfbeta_thresh <- 2 * nrow(dat)^-0.5
		cooksthresh1 <- 4 / (nrow(dat) - 2)
		cooksthresh2 <- 4 / (nrow(dat) - 3)
		inf1 <- influence.measures(ruck$lmod_ivw)$infmat
		inf2 <- influence.measures(ruck$lmod_egger)$infmat
		metrics$dfb1_ivw <- sum(inf1[,1] > dfbeta_thresh) / nrow(dat)
		metrics$dfb2_ivw <- sum(inf1[,2] > dfbeta_thresh) / nrow(dat)
		metrics$dfb3_ivw <- sum(inf1[,3] > dfbeta_thresh) / nrow(dat)
		metrics$cooks_ivw <- sum(inf1[,4] > cooksthresh1) / nrow(dat)
		metrics$dfb1_egger <- sum(inf2[,1] > dfbeta_thresh) / nrow(dat)
		metrics$dfb2_egger <- sum(inf2[,2] > dfbeta_thresh) / nrow(dat)
		metrics$dfb3_egger <- sum(inf2[,3] > dfbeta_thresh) / nrow(dat)
		metrics$cooks_egger <- sum(inf2[,4] > cooksthresh2) / nrow(dat)

		# Homoscedasticity
		metrics$homosc_ivw <- car::ncvTest(ruck$lmod_ivw)$ChiSquare
		metrics$homosc_egg <- car::ncvTest(ruck$lmod_egger)$ChiSquare

		# Normality of residuals
		metrics$shap_ivw <- shapiro.test(residuals(ruck$lmod_ivw))$statistic
		metrics$shap_egger <- shapiro.test(residuals(ruck$lmod_egger))$statistic
		metrics$ks_ivw <- ks.test(residuals(ruck$lmod_ivw), "pnorm")$statistic
		metrics$ks_egger <- ks.test(residuals(ruck$lmod_egger), "pnorm")$statistic

	}
	return(metrics)
}



get_metrics <- function(dat)
{
	metrics <- system_metrics(dat)
	# Steiger
	l0 <- list()
	for(i in 1:nrow(dat))
	{
		l0[[i]] <- steiger_simple(
			dat$pval.exposure[i], 
			dat$pval.outcome[i], 
			dat$samplesize.exposure[i], 
			dat$samplesize.outcome[i]
		)
	}
	steiger_keep <- sapply(l0, function(x) x$dir & x$pval < 0.05)
	metrics$st_correct <- sum(steiger_keep) / nrow(dat)
	metrics$st_unknown <- sum(sapply(l0, function(x) x$pval < 0.05)) / nrow(dat)
	metrics$st_incorrect <- sum(sapply(l0, function(x) !x$dir & x$pval < 0.05)) / nrow(dat)

	dat2 <- dat[steiger_keep, ]
	if(nrow(dat2) > 0 & nrow(dat2) != nrow(dat))
	{
		metrics2 <- system_metrics(dat2)
		names(metrics2) <- paste0(names(metrics2), "_after_steiger")
		metrics <- c(metrics, metrics2)
	}
	return(metrics)
}
