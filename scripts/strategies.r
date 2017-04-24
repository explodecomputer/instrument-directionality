
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

nid <- 500000
eff <- 0.015
g <- rbinom(nid, 2, 0.5)
y <- scale(g) * sqrt(eff) + rnorm(nid, sd=sqrt(1-eff))

cc1 <- y
cc1[y < quantile(y, 0.5)] <- 1
cc1[y >= quantile(y, 0.5)] <- 0

cc2 <- y
cc2[y < quantile(y, 0.1)] <- 1
cc2[y >= quantile(y, 0.1)] <- 0

cc3 <- y
cc3[y < quantile(y, 0.05)] <- 1
cc3[y >= quantile(y, 0.05)] <- 0

cor(g, y)^2
cor(g, cc1)^2
cor(g, cc2)^2
cor(g, cc3)^2

m1 <- glm(cc1 ~ g, family="binomial")
m2 <- glm(cc2 ~ g, family="binomial")
m3 <- glm(cc3 ~ g, family="binomial")

NagelkerkeR2(m1)
NagelkerkeR2(m2)
NagelkerkeR2(m3)


p <- 0.5
1 - (p^p * (1-p)^(1-p))^2
1 - exp(-m1$null/nid)

1 - (p^p * (1-p)^(1-p))^2 - 1
- exp(-m1$null/nid)

log(-(1 - (p^p * (1-p)^(1-p))^2 - 1))
- m1$null/nid


log(-(1 - (p^p * (1-p)^(1-p))^2 - 1)) * -nid
m1$null

or <- log(m1$coefficients[1])

or

f1 <- glm(cc1 ~ sample(g), family="binomial")
f1$dev
f1$null

log(f1$null)

m1
m2
m3
(m1^2 * 2 * 0.5 * 0.5 / var(cc1))
(m2^2 * 2 * 0.5 * 0.5 / var(cc2))
(m3^2 * 2 * 0.5 * 0.5 / var(cc3))

0.5/0.5

dnorm(0.1)^2


log(m1$null)


1-exp(-n1$null/nid)
NagelkerkeR2(m3)

(1 - NagelkerkeR2(m3)$R2) ^ (1/(2/nid))

cor()


index1 <- cc1 == 1
index1[sample(which(cc1 == 0), sum(index1), replace=FALSE)] <- TRUE
table(index1, cc1)

index2 <- cc2 == 1
index2[sample(which(cc2 == 0), sum(index2), replace=FALSE)] <- TRUE
table(index2, cc2)

index3 <- cc3 == 1
index3[sample(which(cc3 == 0), sum(index3), replace=FALSE)] <- TRUE
table(index3, cc3)



(m1c <- glm(cc1[index1] ~ g[index1], family="binomial"))
(m2c <- glm(cc2[index2] ~ g[index2], family="binomial")$coefficients[2])
(m3c <- glm(cc3[index3] ~ g[index3], family="binomial")$coefficients[2])

sum(g[index1]) / (2 * sum(index1))
sum(g[index2]) / (2 * sum(index2))
sum(g[index3]) / (2 * sum(index3))

(m1c^2 * 2 * 0.5 * 0.5 / var(cc1))^2
cor()


NagelkerkeR2()




set.seed(17)
i<-0
sim <- replicate(1e4, {
  while(TRUE) {
    x <- matrix(round(rexp(4), 2), 2, 2)
    if(all(rowSums(x) > 0) && all(colSums(x) > 0) && x[1,2]*x[2,1] > 0) break
  }
  x <- x / sum(x)
  beta <- rowSums(x)[1]
  gamma <- colSums(x)[1]
  rho <- x[1,1]*x[2,2] / (x[1,2]*x[2,1])

  y<- f(beta, gamma, rho)
  delta <- try(zapsmall(c(1, sqrt(crossprod(as.vector(x-y)))))[2])
  if ("try-error" %in% class(delta)) cat("Error processing ", x, "\n")
  delta
})
max(sim)

f(1.2, 0.2, 0.5)
f(log(1.2), 0.2, 0.5)

A <- 200
B <- 400
C <- 200
D <- 1000
OR <- A*D/(B*C)
N <- A+B+C+D
f((A+B)/N, (A+C)/N, OR)




A/N
B/N
C/N
D/N


cont_from_or <- function(prevalence, allele_frequency, oddsratio, eps=1e-15)
{
	a <- oddsratio-1
	b <- (prevalence+allele_frequency)*(1-oddsratio)-1
	c_ <- oddsratio*prevalence*allele_frequency

	if (abs(a) < eps) {
	z <- -c_ / b
	} else {
	d <- b^2 - 4*a*c_
	if (d < eps*eps) s <- 0 else s <- c(-1,1)
	z <- (-b + s*sqrt(max(0, d))) / (2*a)
	}
	y <- vapply(z, function(a) zapsmall(matrix(c(a, allele_frequency-a, prevalence-a, 1+a-prevalence-allele_frequency), 2, 2)), matrix(0.0, 2, 2))
	i <- apply(y, 3, function(u) all(u >= 0))
	return(y[,,i])
}

cont_from_or(0.5, 0.5, exp(m1$coefficients[2]))
sum(cc1 & )/nid

(sum(cc1 & g==2)*2 + sum(cc1 & g==1)) / (2*nid)


