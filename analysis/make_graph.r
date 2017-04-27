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

m_mode <- lapply(d, function(x) {
	group_by(x, id.exposure, id.outcome) %>%
	do({
		x <- .
		print(head(x))
		mr_mode(x)
	})
})




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


