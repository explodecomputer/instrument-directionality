library(tidyverse)
library(TwoSampleMR)

load("../data/extract_everything2.rdata")
load("../data/extract_data.rdata")

# Identify NULL elements in list

m_ivw <- lapply(d, function(x) {
	mr(x, method_list=c("mr_wald_ratio", "mr_ivw"))
})

m_weighted_median <- lapply(d, function(x) {
	mr(x, method_list=c("mr_wald_ratio", "mr_weighted_median"))
})

m_all <- lapply(d, function(x) {
	x <- subset(x, id.exposure != id.outcome)
	group_by(x, exposure, outcome, id.exposure, id.outcome) %>%
	do(
	{
		if(nrow(x) > 3)
		{
			res <- try(mr_all(.)$res)
			if(class(res) != "try-error")
			{
				return(res)
			} else {
				return(NULL)
			}
		}
	})
})

mr_all(x[1:14,])


a <- extract_instruments(300)
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a,b)

mr_mode(dat)

# Mode
# Median
# IVW
# Egger
# Rucker
# Rucker (CD)

# Steiger




make_tab <- function(m)
{
	index <- which(sapply(m, function(x) nrow(x)==0))
	index <- which(sapply(m, function(x) nrow(x)==0))
	temp <- m[-index]
	m <- bind_rows(temp)
	# m <- subset(m, id.exposure != id.outcome)
	m$pval_bonf <- p.adjust(m$pval, "bonferroni")
	m$pval_fdr <- p.adjust(m$pval, "fdr")
	return(m)
}

mwm <- make_tab(m_weighted_median)
mivw <- make_tab(m_ivw)


table(mwm$pval < 0.05)
table(mivw$pval < 0.05)

table(mwm$pval_fdr < 0.05)
table(mivw$pval_fdr < 0.05)

table(mwm$pval_bonf < 0.05)
table(mivw$pval_bonf < 0.05)

tail(subset(mwm, pval_bonf < 0.05))


for(i in 1:nrow(traits2))
{
	temp[[i]] <- subset(temp[[i]], !id.outcome == traits2$id[i])
}

m_ivwd <- bind_rows(temp)
dim(m_ivwd)
head(m_ivwd)
(m_ivwd$pval < 0.05/nrow(m_ivwd), na.rm=TRUE) 

a <- subset(m_ivwd, pval < (0.05/nrow(m_ivwd)))
m_ivwd$fdr <- p.adjust(m_ivwd$pval, "fdr")

head(m_ivwd)
a <- subset(m_ivwd, fdr < 0.05)

dim(a)
table(a$exposure)

subset(a, exposure == "Intracranial volume || ENIGMA || 2015 || mm3")



steiger_simple <- function(p_exp, p_out, n_exp, n_out)
{
	r_exp <- sqrt(sum(get_r_from_pn(p_exp, n_exp)^2))
	r_out <- sqrt(sum(get_r_from_pn(p_out, n_out)^2))
	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	return(list(dir=r_exp > r_out, pval=rtest$p))
}


steiger_dat <- function(dat, steiger_threshold)
{
	group_by(dat, id.exposure, id.outcome)
}


sapply(m, function(x) length(is.na(x$samplesize.outcome)))






###


library(plyr)
load("../data/m_all.rdata")
load("../data/mr_res.rdata")
load("../data/extract_data.rdata")
nom <- lapply(d, function(x1)
{
	x1 <- subset(x1, id.exposure != id.outcome)
	dlply(x1, .(exposure, outcome, id.exposure, id.outcome), function(x)
	{
		if(nrow(x) > 3) return(x)
		})
	})
stopifnot(all(sapply( m_all, length) == sapply(nom, length)))
l <- list()
k <- 1
for(i in 1:length(m_all))
{
	message(i)
	a <- attributes(m_all[[i]])$split_labels
	if(!is.null(a))
	{
		b <- attributes(m_all[[i]])$names
		a$lab <- with(a, paste(id.exposure, id.outcome, exposure, outcome, sep="."))
		index <- match(a$lab, b)
		a <- a[index, ]
		stopifnot(all(a$lab == b))
		for(j in 1:nrow(a))
		{
			if(j %in% index)
			{
				m_all[[i]][[j]]$exposure <- a$exposure[j]
				m_all[[i]][[j]]$outcome <- a$outcome[j]
				m_all[[i]][[j]]$id.exposure <- a$id.exposure[j]
				m_all[[i]][[j]]$id.outcome <- a$id.outcome[j]
				if(is.data.frame(m_all[[i]][[j]])) {
					if(nrow(m_all[[i]][[j]]) == 17){
						l[[k]] <- m_all[[i]][[j]]
						k <- k+1
					}
				}
			}
		}
	}
}

length(l)
table(sapply(l, ncol))
m_plus3 <- bind_rows(l)
mivwa <- subset(mivw, nsnp <= 3)
names(mivwa) <- c("id.exposure", "id.outcome", "outcome", "exposure", "Method", "nsnp", "Estimate", "SE", "P", "pval_bonf", "pval_fdr")
mivwa$CI_low <- mivwa$Estimate - 1.96 * mivwa$SE
mivwa$CI_upp <- mivwa$Estimate + 1.96 * mivwa$SE

m <- unique(m_plus3$Method)
m_plus3$pval_bonf <- NA
m_plus3$pval_fdr <- NA
for(i in 1:length(m))
{
	index <- m_plus3$Method == m[i]
	pval <- c(m_plus3$P[index], mivwa$P)
	m_plus3$pval_bonf[index] <- p.adjust(pval, method="bonferroni")[1:sum(index)]
	m_plus3$pval_fdr[index] <- p.adjust(pval, method="fdr")[1:sum(index)]
}


save(m_plus3, mivwa, file="../data/mres_organised.rdata")


steiger_simple <- function(p_exp, p_out, n_exp, n_out)
{
	p_exp[p_exp == 0] <- 1e-200
	p_out[p_out == 0] <- 1e-200
	r_exp <- sqrt(sum(TwoSampleMR::get_r_from_pn(p_exp, n_exp)^2))
	r_out <- sqrt(sum(TwoSampleMR::get_r_from_pn(p_out, n_out)^2))
	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	return(list(dir=r_exp > r_out, pval=rtest$p))
}



load("../data/mres_organised.rdata")
load("../data/extract_data.rdata")





instrument_directionality <- function(dat, thresh=0.05)
{
	# Choose only traits that are continuous
	# Need pval and sample size
	message(dat$exposure[1], ": ", nrow(dat), " associations")
	dat <- subset(dat, units.exposure != "log odds" & units.outcome != "log odds")
	dat <- subset(dat, !is.na(samplesize.exposure) & !is.na(pval.exposure) & !is.na(samplesize.outcome) & !is.na(pval.outcome))
	message(dat$exposure[1], ": ", nrow(dat), " associations left")
	if(nrow(dat) > 0)
	{
		st <- list()
		for(i in 1:nrow(dat))
		{
			st[[i]] <- steiger_simple(dat$pval.exposure[i], dat$pval.outcome[i], dat$samplesize.exposure[i], dat$samplesize.outcome[i])
		}
		dat$dir <- sapply(st, function(x) x$dir)
		dat$st_pval <- sapply(st, function(x) x$pval)
		dat <- subset(dat, dir & st_pval < thresh)
		message(dat$exposure[1], ": ", nrow(dat), " positive instruments left")
		return(dat)
	} else {
		message("Insufficient data")
		return(NULL)
	}
}

ds <- lapply(d, instrument_directionality)
dat <- bind_rows(ds)
traits <- subset(traits, trait %in% c(dat$exposure, dat$outcome))
save(dat, traits, file="../data/continuous_steiger_dat.rdata")



##

library(TwoSampleMR)
library(dplyr)
load("../data/continuous_steiger_dat.rdata")

dat <- subset(dat, id.exposure != id.outcome)
a <- mr_all(dat)
b <- group_by(a, Method) %>%
mutate(pval_bonf = p.adjust(P, method="bonferroni"), pval_fdr = p.adjust(P, method="fdr"))

group_by(b, Method) %>%
summarise(bonf=sum(pval_bonf < 0.05, na.rm=TRUE), fdr = sum(pval_fdr < 0.05, na.rm=TRUE))
m_all_st <- b
save(m_all_st, file="../data/m_all_st.rdata")


b$lab <- paste0(b$exposure, b$outcome, sep=" - ")


##


source("~/repo/graph_mr/scripts/functions.r")
load("../data/mres_organised.rdata")
mivwa <- subset(mivwa, !is.na(P))


extract_result <- function(m_plus3, mivwa, method)
{
	a <- subset(m_plus3, Method == method)
	return(bind_rows(a, mivwa))
}

get_all_paths <- function(mres, thresh, from)
{
	g1 <- filter(mres, P < thresh) %>%
		select(exposure, outcome) %>%
		as.matrix %>%
		graph_from_edgelist(directed=TRUE)
	res <- all_shortest_paths(g1, from)
	return(res)	
}


make_matrix <- function(mres)
{
	traits <- unique(c(mres$exposure, mres$outcome))
	mat <- matrix(NA, length(traits), length(traits))
	i <- cbind(
		match(mres$exposure, traits),
		match(mres$outcome, traits)
	)
	nas <- apply(i, 1, function(x) any(is.na(x)))
	mat[i] <- mres$P
	colnames(mat) <- traits
	rownames(mat) <- traits

	print(mat[traits == "LDL cholesterol", traits == "Coronary heart disease"])
	print(subset(mres, exposure == "LDL cholesterol" & outcome == "Coronary heart disease")$P)

	stopifnot(mat[1, 100] == subset(mres, exposure == traits[1] & outcome == traits[100])$P)
	mat <- t(mat)
	print(mat[traits == "Coronary heart disease", traits == "LDL cholesterol"])
	print(subset(mres, exposure == "LDL cholesterol" & outcome == "Coronary heart disease")$P)
	return(t(mat))
}

find_path <- function(mat, thresh, from, to, minn=2, maxn=Inf)
{
	g1 <- graph_from_adjacency_matrix(mat < thresh, mode="directed")
	res <- all_shortest_paths(g1, from, to)$res
	
	a <- length(res) > 0
	if(a)
	{
		res <- res[sapply(res, function(x) length(x) >= minn & length(x) <= maxn)]
		if(length(res) == 0) return(1)
		n <- length(res[[1]]) - 1
		c1 <- match(names(res[[1]])[1:n], colnames(mat))
		c2 <- match(names(res[[1]])[2:(n+1)], colnames(mat))
		return(max(mat[cbind(c1, c2)]))
	} else {
		return(1)
	}
}

mat <- make_matrix(mwmed)
g1 <- graph_from_adjacency_matrix(mat < 1e-4, mode="directed")

find_path(mat, 1e-4, from="Body mass index", to="Coronary heart disease")
find_path(mat, 1e-4, "LDL cholesterol",  "Coronary heart disease")
find_path(mat, 1e-4, "Coronary heart disease", "LDL cholesterol")


matp <- permute_matrix(mat)
find_path(matp, 1e-4, "LDL cholesterol",  "Coronary heart disease")


get_empirical_pval <- function(mat, thresh, from, to, nboot=1000, minn=3, maxn=Inf)
{
	p <- find_path(mat, thresh, from, to)
	out <- rep(0, nboot)
	for(i in 1:nboot)
	{
		out[i] <- find_path(permute_matrix(mat), thresh, from, to, minn, maxn)
	}
	return(list(
		empirical = sum(p > out) / nboot,
		original = p,
		permutations = out
	))
}


a <- get_empirical_pval(mat, fdrthresh, "Body mass index",  "Coronary heart disease")


# How many FDR < 0.05 and bonferroni < 0.05 using each different method

group_by(m_plus3, Method) %>%
summarise(bonf = sum(pval_bonf < 0.05) + sum(mivwa$pval_bonf < 0.05), fdr = sum(pval_fdr < 0.05) + sum(mivwa$pval_fdr < 0.05))

sum(mivwa$pval_bonf < 0.05)
sum(mivwa$pval_fdr < 0.05)




# What thresholds to use for different path lengths
fdrthresh <- 0.0017
mwmed <- extract_result(m_plus3, mivwa, "Weighted median")
mwmode <- extract_result(m_plus3, mivwa, "Weighted mode")
calc_number_of_paths(150, 0.01, 5) %>% colSums
calc_number_of_paths(150, 0.01) %>% colSums
calc_number_of_paths(150, 0.0017) %>% colSums

# Seems like a link threshold of FDR (0.05) = p(0.0017) is fine for preventing big ramp up with 150 traits

get_all_graph_outcomes <- function(mres, thresh, from)
{
	res <- get_all_paths(mres, thresh, from)
	res <- res$res[sapply(res$res, function(x) !is.null(x))]
	print(length(res))
	dat <- tibble(
		exposure = from,
		outcome = sapply(res, function(x) last(names(x))),
		links = sapply(res, length),
		path = sapply(res, function(x) paste(names(x), collapse=" -> "))
	)
	dat <- dat[order(dat$links), ]
	dat <- subset(dat, !duplicated(outcome))
	dat <- subset(dat, links != 1)
	return(dat)
}

a <- get_all_paths(mwmed, fdrthresh, "LDL cholesterol")

l <- list()
exposures <- unique(subset(mwmed, pval_fdr < 0.05)$exposure)
for(i in 1:length(exposures))
{
	l[[i]] <- get_all_graph_outcomes(mwmed, fdrthresh, exposures[i])
	print(dim(l[[i]]))
}
l <- bind_rows(l)





group_by(mwmode, exposure) %>%
summarise(fdr = sum(pval_fdr < 0.05), bonf = sum(pval_bonf < 0.05))

plot(-log10(mwmed$P), -log10(mwmed$pval_fdr), xlim=c(0,5), ylim=c(0,5))
plot(-log10(mwmed$P), -log10(mwmed$pval_bonf), xlim=c(0,5), ylim=c(0,5))



mat <- find_path(mivw, 0.05/16000, "LDL cholesterol", "Coronary heart disease")


find_path <- function(mres, thresh, from, to)
{
	g1 <- filter(mres, P < thresh) %>%
		select(exposure, outcome) %>%
		as.matrix %>%
		graph_from_edgelist(directed=TRUE)
	res <- all_shortest_paths(g1, from)
	a <- length(res$res) > 0
	if(a)
	{
		n <- length(res$res[[1]]) - 1
		c1 <- res$res[[1]][1:n]
		c2 <- res$res[[1]][2:(n+1)]
		return(max(mat[rbind(c1, c2)]))
	} else {
		return(1)
	}
}




a <- get_all_paths(mivw, 1e-4, "Years of schooling")
a <- get_all_paths(mwm, 1e-4, "Intracranial volume")


find


g <- make_ring(10)
     distances(g)
     shortest_paths(g, 5)
     all_shortest_paths(g, 1, 6:8)

library(plyr)
m_all <- lapply(d, function(x1) {
	x1 <- subset(x1, id.exposure != id.outcome)
	plyr::dlply(x1, .(id.exposure, id.outcome, exposure, outcome), function(x) {
		if(nrow(x) > 3)
		{
			res <- try(mr_all(x)$res)
			if(class(res) != "try-error")
			{
				return(res)
			} else {
				return(NULL)
			}

		} else {
			return(NULL)
		}
	})
})



