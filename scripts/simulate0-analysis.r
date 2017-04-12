library(tidyverse)

load("../results/simulate0.rdata")

dat_sel$method <- "Rucker"
dat_sel$var_g1y_lab <- paste0("sigma(a) = ", dat_sel$var_g1y)
dat_sel$mu_g1y_lab <- paste0("mu(a) = ", dat_sel$mu_g1y)

dats <- rbind(as.data.frame(dat_sel), dat)
dats$var_g1y_lab <- paste0("sigma(a) = ", dats$var_g1y)
dats$mu_g1y_lab <- paste0("mu(a) = ", dats$mu_g1y)


# Bias not null
ggplot(filter(dats, var_xy != 0), aes(x=as.factor(nsnp.x), y=var_xy-b)) +
geom_boxplot(aes(fill=method)) +
facet_grid(var_g1y ~ mu_g1y) +
scale_fill_brewer(type="qual") +
labs(x="Number of instruments", y="Bias")

# Rucker choice not null
ggplot(filter(dat_sel, var_xy != 0), aes(x=as.factor(nsnp.x))) +
geom_bar(aes(fill=rucker)) +
facet_grid(var_g1y_lab ~ mu_g1y_lab) +
scale_fill_brewer(type="qual") +
labs(x="Number of instruments", y="Frequency")

# Rucker choice under null
ggplot(filter(dat_sel, var_xy==0), aes(x=as.factor(nsnp.x))) +
geom_bar(aes(fill=rucker)) +
facet_grid(var_g1y_lab ~ mu_g1y_lab) +
scale_fill_brewer(type="qual") +
labs(x="Number of instruments", y="Frequency")

# FDR

dats_fdr <- filter(dats, var_xy==0) %>%
	group_by(nsnp.x, method, var_g1y, mu_g1y) %>%
	summarise(fdr = sum(pval < 0.05) / n())

ggplot(dats_fdr, aes(x=method, y=fdr)) +
geom_point(aes(colour=as.factor(nsnp.x)), size=2) +
facet_grid(var_g1y ~ mu_g1y) +
scale_colour_brewer() +
labs(x="Number of instruments", y="FDR") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


# Power
dats_tdr <- filter(dats, var_xy!=0) %>%
	group_by(nsnp.x, method, var_g1y, mu_g1y, var_xy) %>%
	summarise(tdr = sum(pval < 0.05) / n())

ggplot(filter(dats_tdr, var_xy==0.1), aes(x=as.factor(nsnp.x), y=tdr)) +
geom_jitter(aes(colour=method, size=var_xy)) +
facet_grid(var_g1y ~ mu_g1y) +
scale_colour_brewer(type="qual") +
labs(x="Number of instruments", y="TDR")

