
steiger_simple <- function(p_exp, p_out, n_exp, n_out)
{
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

