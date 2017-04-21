library(dplyr)

makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	sc <- sum(cors^2)
	if(sc >= 1)
	{
		print(sc)
		stop("effects explain more than 100% of variance")
	}
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

chooseEffects <- function(nsnp, totvar, sqrt=TRUE, mua=0)
{
	eff <- rnorm(nsnp)
	eff <- sign(eff) * eff^2
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out + mua)
}

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, n=n
	))
}

gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 6)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

get_effs <- function(x, y, g)
{
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)

	dat <- data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=gwasx$bhat,
		beta.outcome=gwasy$bhat,
		se.exposure=gwasx$se,
		se.outcome=gwasy$se,
		pval.exposure=gwasx$pval,
		pval.outcome=gwasy$pval,
		samplesize.exposure=gwasx$n,
		samplesize.outcome=gwasy$n,
		mr_keep=TRUE
	)
	return(dat)
}

dat_sign <- function(dat)
{
	sign0 <- function(x) {
		x[x == 0] <- 1
		return(sign(x))
	}
	index <- sign0(dat$beta.exposure) == -1
	dat$beta.exposure <- abs(dat$beta.exposure)
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
	return(dat)
}

recode_dat <- function(dat)
{
	a <- lm(beta.outcome ~ beta.exposure, dat)$coefficients[1]
	index <- dat$beta.exposure < 0
	dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1 + 2 * a
	dat$index <- index
	return(dat)
}

make_effs <- function(ninst1, var_g1u=0, var_g1x, var_g1y=0, mu_g1y=0, var_xy, var_ux=0, var_uy=0, ninst2=0, var_g2u=0, var_g2x=0, var_g2y=0, mu_g2x=0, ninstu=0, var_guu=0)
{
	# 1 SNP influences one confounder
	eff_g1u <- chooseEffects(ninst1, var_g1u)
	eff_g1x <- chooseEffects(ninst1, var_g1x)
	eff_g1y <- chooseEffects(ninst1, var_g1y, mua=mu_g1y)

	eff_g2u <- chooseEffects(ninst2, var_g2u)
	eff_g2x <- chooseEffects(ninst2, var_g2x, mua=mu_g2x)
	eff_g2y <- chooseEffects(ninst2, var_g2y)

	eff_guu <- chooseEffects(ninstu, var_guu)

	eff_xy <- var_xy
	eff_ux <- var_ux
	eff_uy <- var_uy

	return(list(
		eff_g1u = eff_g1u,
		eff_g1x = eff_g1x,
		eff_g1y = eff_g1y,
		eff_g2u = eff_g2u,
		eff_g2x = eff_g2x,
		eff_g2y = eff_g2y,
		eff_ux = eff_ux,
		eff_uy = eff_uy,
		eff_xy = eff_xy,
		eff_guu = eff_guu
	))
}

make_pop <- function(effs, nid)
{
	ninst1 <- length(effs$eff_g1x)
	ninst2 <- length(effs$eff_g2x)
	ninstu <- length(effs$eff_guu)
	G1 <- matrix(rbinom(nid * ninst1, 2, 0.5), nid, ninst1)
	G2 <- matrix(rbinom(nid * ninst2, 2, 0.5), nid, ninst2)
	Gu <- matrix(rbinom(nid * ninstu, 2, 0.5), nid, ninstu)

	u <- makePhen(c(effs$eff_g1u, effs$eff_g2u, effs$eff_guu), cbind(G1, G2, Gu))
	x <- makePhen(c(effs$eff_g1x, effs$eff_g2x, effs$eff_ux), cbind(G1, G2, u))
	y <- makePhen(c(effs$eff_xy, effs$eff_uy, effs$eff_g1y, effs$eff_g2y), cbind(x, u, G1, G2))
	return(list(
		x=x,
		y=y,
		u=u,
		G1=G1,
		G2=G2,
		Gu=Gu
	))
}

make_dat <- function(exposure, outcome)
{
	dat <- cbind(
		exposure[,grepl("exposure", names(exposure))],
		outcome[,grepl("outcome", names(outcome))]
	)
	dat$mr_keep <- TRUE
	return(dat)
}

analyse_simulation <- function(dat, pop2)
{
	prs <- pop2$G %*% sign(dat$beta.exposure)
	wprs <- pop2$G %*% dat$beta.exposure
	mod1 <- coefficients(summary(lm(pop2$y ~ prs)))
	mod2 <- coefficients(summary(lm(pop2$y ~ wprs)))
	scoreres <- data.frame(
		id.exposure = "X", id.outcome = "Y", exposure = "X", outcome = "Y",
		method = c("Basic PRS", "Weighted PRS"),
		nsnp = nrow(dat),
		b = c(mod1[2,1], mod2[2,1]),
		se = c(mod1[2,2], mod2[2,2]),
		pval = c(mod1[2,4], mod2[2,4])
	)
	ivres <- mr(dat, method_list=mr_method_list()$obj[c(10, 6, 8)])
	res <- rbind(ivres, scoreres)
	return(res)
}

make_datg <- function(gwasx, gwasy)
{
	dat <- data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=gwasx$bhat,
		beta.outcome=gwasy$bhat,
		se.exposure=gwasx$se,
		se.outcome=gwasy$se,
		pval.exposure=gwasx$pval,
		pval.outcome=gwasy$pval,
		samplesize.exposure=gwasx$n,
		samplesize.outcome=gwasy$n,
		mr_keep=TRUE,
		inst=gwasx$inst
	)
	return(dat)

}

init_parameters <- function(nsnp_x, var_gx.x, var_x.y, var_gx.y=0, nsnp_y=0, var_gy.y=0, mu_gx.y=0, var_gy.x=0, mu_gy.x=0)
{
	parameters <- list(
		# Causal effect
		var_x.y = var_x.y,

		# Direct effects on x
		nsnp_x = nsnp_x,
		var_gx.x = var_gx.x,
		var_gx.y = var_gx.y,
		mu_gx.y = mu_gx.y,

		nsnp_y = nsnp_y,
		var_gy.y = var_gy.y,
		var_gy.x = var_gy.x,
		mu_gy.x = mu_gy.x,
		u = list()
	)
	return(parameters)
}

add_u <- function(parameters, nsnp_u, var_u.x, var_u.y, var_gu.u)
{
	nom <- "u"
	if(var_u.x == 0 & var_u.y != 0) nom <- "b"
	if(var_u.x != 0 & var_u.y == 0) nom <- "a"
	i <- length(parameters$u) + 1
	parameters$u[[i]] <- list(
		nsnp_u = nsnp_u,
		var_u.x = var_u.x,
		var_u.y = var_u.y,
		var_gu.u = var_gu.u,
		name_u = nom
	)
	return(parameters)
}

generate_system_effects <- function(parameters)
{
	nu <- length(parameters$u)
	if(nu > 0)
	{
		for(i in 1:nu)
		{
			parameters$u[[i]]$eff_gu.u <- chooseEffects(parameters$u[[i]]$nsnp_u, parameters$u[[i]]$var_gu.u)
			parameters$u[[i]]$eff_u.x <- chooseEffects(1, parameters$u[[i]]$var_u.x)
			parameters$u[[i]]$eff_u.y <- chooseEffects(1, parameters$u[[i]]$var_u.y)
		}
	}

	parameters$eff_gx.x <- chooseEffects(parameters$nsnp_x, parameters$var_gx.x)
	parameters$eff_gx.y <- chooseEffects(parameters$nsnp_x, parameters$var_gx.y, mua=parameters$mu_gx.y)

	parameters$eff_gy.x <- chooseEffects(parameters$nsnp_y, parameters$var_gy.x, mua=parameters$mu_gy.x)

	parameters$eff_gy.y <- chooseEffects(parameters$nsnp_y, parameters$var_gy.y)

	parameters$eff_u.x <- sapply(parameters$u, function(x) chooseEffects(1, x$var_u.x))

	parameters$eff_u.y <- sapply(parameters$u, function(x) chooseEffects(1, x$var_u.y))
	parameters$eff_x.y <- sqrt(parameters$var_x.y)

	return(parameters)
}

simulate_population <- function(parameters, nid)
{
	require(dplyr)

	Gx <- matrix(rbinom(parameters$nsnp_x * nid, 2, 0.5), nid, parameters$nsnp_x)

	Gy <- matrix(rbinom(parameters$nsnp_y * nid, 2, 0.5), nid, parameters$nsnp_y)

	U <- lapply(parameters$u, function(param)
	{
		G <- matrix(rbinom(param$nsnp_u * nid, 2, 0.5), nid, param$nsnp_u)
		u <- makePhen(param$eff_gu.u, G)
		return(list(p=u, G=G, nom=param$name_u))
	})


	bx <- parameters$eff_gx.x
	by <- parameters$eff_gx.y

	Gt <- Gx
	if(parameters$nsnp_y > 0)
	{
		bx <- c(bx, parameters$eff_gy.x)
		by <- c(by, parameters$eff_gy.y)
		Gt <- cbind(Gt, Gy)
	}

	if(length(parameters$u) > 0)
	{
		bx <- c(bx, sapply(parameters$u, function(x) x$eff_u.x))
		by <- c(by, sapply(parameters$u, function(x) x$eff_u.y))
		Gt <- cbind(Gt, do.call(cbind, lapply(U, function(x) x$p)))
	}
	by <- c(by, parameters$eff_x.y)

	x <- makePhen(bx, Gt)
	y <- makePhen(by, cbind(Gt, x))

	return(list(
		y=y,
		x=x,
		Gx=Gx,
		Gy=Gy,
		U=U
	))

}


system_effs <- function(sim)
{
	# Get effects of all X SNPs

	l <- list()

	gx.x <- gwas(sim$x, sim$Gx)
	gx.x$inst <- "x"
	gx.x$snp <- 1:nrow(gx.x)
	
	gx.y <- gwas(sim$y, sim$Gx)
	gx.y$inst <- "x"
	gx.y$snp <- 1:nrow(gx.y)

	gx <- gx.x
	gy <- gx.y

	# Get effects of all SNPs on Y

	if(ncol(sim$Gy) > 0)
	{
		gy.x <- gwas(sim$x, sim$Gy)
		gy.x$inst <- "y"
		gy.x$snp <- 1:nrow(gy.x)
		gx <- rbind(gx, gy.x)

		gy.y <- gwas(sim$y, sim$Gy)
		gy.y$inst <- "y"
		gy.y$snp <- 1:nrow(gy.y)
		gy <- rbind(gy, gy.y)
	}

	nconf <- length(sim$U)
	if(nconf > 0)
	{
		gu.x <- list()
		gu.y <- list()
		gx.u <- list()
		gy.u <- list()
		gu.u <- list()
		for(i in 1:nconf)
		{
			gu.u[[i]] <- gwas(sim$U[[i]]$p, sim$U[[i]]$G)
			gu.u[[i]]$inst <- paste0("u", i)
			gu.u[[i]]$snp <- 1:nrow(gu.u[[i]])

			gu.x[[i]] <- gwas(sim$x, sim$U[[i]]$G)
			gu.x[[i]]$inst <- paste0("u", i)
			gu.x[[i]]$snp <- 1:nrow(gu.x[[i]])

			gx.u[[i]] <- gwas(sim$U[[i]]$p, sim$Gx)
			gx.u[[i]]$inst <- "x"
			gx.u[[i]]$snp <- 1:nrow(gx.u[[i]])

			gu.y[[i]] <- gwas(sim$y, sim$U[[i]]$G)
			gu.y[[i]]$inst <- paste0("u", i)
			gu.y[[i]]$snp <- 1:nrow(gu.y[[i]])

			if(ncol(sim$Gy) > 0)
			{
				gy.u[[i]] <- gwas(sim$U[[i]]$p, sim$Gy)
				gy.u[[i]]$inst <- "y"
				gy.u[[i]]$snp <- 1:nrow(gy.u[[i]])
			}
		}
		gx <- rbind(gx, bind_rows(gu.x))
		gy <- rbind(gy, bind_rows(gu.y))
		gu <- list()
		for(i in 1:nconf)
		{
			if(ncol(sim$Gy) > 0)
			{
				gu[[i]] <- rbind(gx.u[[i]], gy.u[[i]], gu.u[[i]])
			} else {
				gu[[i]] <- rbind(gx.u[[i]], gy.u[[i]])
			}
		}
		names(gu) <- paste0("u", 1:nconf)
		l <- gu
	}
	l$x <- gx
	l$y <- gy
	return(l)
}

system_dat <- function(gwasx, gwasy)
{
	gwasx <- filter(gwasx, inst %in% gwasy$inst)
	gwasy <- filter(gwasy, inst %in% gwasx$inst)
	dat <- data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=gwasx$bhat,
		beta.outcome=gwasy$bhat,
		se.exposure=gwasx$se,
		se.outcome=gwasy$se,
		pval.exposure=gwasx$pval,
		pval.outcome=gwasy$pval,
		samplesize.exposure=gwasx$n,
		samplesize.outcome=gwasy$n,
		mr_keep=TRUE,
		inst=gwasx$inst
	)
	return(recode_dat(dat))
}

create_system <- function(nidx, nidy, nidu, nu, na, nb, var_x.y, var_gx.x=0.02, var_gy.y=0.02)
{
	parameters <- init_parameters(nsnp_x=20, nsnp_y=20, var_gx.x=var_gx.x, var_gy.y=var_gy.y, var_x.y=var_x.y)
	if(nu > 0)
	{
		for(i in 1:nu)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=runif(1, min=0.01, max=0.1), 
				var_u.y=runif(1, min=0.01, max=0.1), 
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}
	if(na > 0)
	{
		for(i in 1:nu)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=runif(1, min=0.01, max=0.1), 
				var_u.y=0,
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}
	if(nb > 0)
	{
		for(i in 1:nu)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=0,
				var_u.y=runif(1, min=0.01, max=0.1),
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}

	parameters <- generate_system_effects(parameters)

	message("X")
	pop <- simulate_population(parameters, nidx)
	x <- system_effs(pop)

	message("Y")
	pop <- simulate_population(parameters, nidy)
	y <- system_effs(pop)

	u <- list()
	if(nu > 0)
	{
		for(i in 1:nu)
		{
			message("U: ", i, " of ", nu)
			pop <- simulate_population(parameters, nidu)
			u[[i]] <- system_effs(pop)
		}
	}

	a <- list()
	if(na > 0)
	{
		for(i in 1:na)
		{
			message("A: ", i, " of ", na)
			pop <- simulate_population(parameters, nidu)
			a[[i]] <- system_effs(pop)
		}
	}

	b <- list()
	if(nb > 0)
	{
		for(i in 1:nb)
		{
			message("B: ", i, " of ", nb)
			pop <- simulate_population(parameters, nidu)
			b[[i]] <- system_effs(pop)
		}
	}
	u <- c(u, a, b)
	names(u) <- paste0("u", 1:length(u))
	return(list(x=x, y=y, u=u, parameters=parameters))
}
