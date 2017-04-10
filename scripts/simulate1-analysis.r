library(tidyverse)
load("../results/simulate1.rdata")

ggplot(res, aes(x=inst, y=-log10(pval))) +
geom_boxplot(aes(fill=method)) +
facet_grid(dir ~ ., scale="free")

filter(res) %>%
group_by(method, dir,inst) %>%
summarise(psig = sum(pval < 0.05)/n())





summary(lm(b ~ var_xy + var_ux + var_uy + var_g1x + var_g2x + var_g1y + var_g2y + var_guu, data=filter(res, method == "mr_egger_regression", inst == "true", dir == "yx")))

summary(lm(b ~ var_xy + var_ux + var_uy + var_g1x + var_g2x + var_g1y + var_g2y + var_guu, data=filter(res, method == "mr_ivw", inst == "true", dir == "yx")))



filter(selection) %>%
	group_by(dir, var) %>%
	summarise(
		prop = mean(included / (included + excluded))
	)

filter(res, method == "mr_egger_regression", inst == "true", dir == "yx") %>% head


resw <- select(res, sim, inst, b, dir, method) %>% spread(key=inst, value=b)

ggplot(filter(resw, method == "mr_ivw"), aes(x=true, y=sel)) +
geom_point() +
facet_grid(dir ~ .)

ggplot(filter(resw, method == "mr_ivw"), aes(x=true, y=all)) +
geom_point() +
facet_grid(dir ~ .)



source("functions.r")
source("simulation-functions.r")


rxy <- array(0, 100)
ryx <- array(0, 100)
for(i in 1:100)
{
	message(i)
	effs <- make_effs(ninst1=50, ninst2=50, var_xy=0.5, var_g1x=0.2, var_g2y=0.2)
	pop1 <- make_pop(effs, 10000)
	pop2 <- make_pop(effs, 10000)
	d1 <- get_effs(pop1$x, pop1$y, cbind(pop1$G1))
	d2 <- get_effs(pop2$x, pop2$y, cbind(pop2$G1))
	dxy <- recode_dat(make_dat(d1, d2)) %>% filter(pval.exposure < 5e-8)
	d1 <- get_effs(pop1$y, pop1$x, cbind(pop1$G2))
	d2 <- get_effs(pop2$y, pop2$x, cbind(pop2$G2))
	dyx <- recode_dat(make_dat(d1, d2)) %>% filter(pval.exposure < 5e-8)
	rxy[i] <- mr(dxy, meth="mr_ivw")$pval
	ryx[i] <- mr(dyx, meth="mr_ivw")$pval
}

rxy_p <- array(0, 100)
ryx_p <- array(0, 100)

p <- 1
while(p > 1e-4)
{
	message(i)
	effs <- make_effs(ninst1=50, ninst2=50, var_xy=0.5, var_g1x=0.2, var_g2y=0.2, var_g1y=0.1, var_g2x=0.1)
	pop1 <- make_pop(effs, 10000)
	pop2 <- make_pop(effs, 10000)
	d1 <- get_effs(pop1$x, pop1$y, cbind(pop1$G1))
	d2 <- get_effs(pop2$x, pop2$y, cbind(pop2$G1))
	dxy <- recode_dat(make_dat(d1, d2)) %>% filter(pval.exposure < 5e-8)
	d1 <- get_effs(pop1$y, pop1$x, cbind(pop1$G2))
	d2 <- get_effs(pop2$y, pop2$x, cbind(pop2$G2))
	dyx <- recode_dat(make_dat(d1, d2)) %>% filter(pval.exposure < 5e-8)
	p <- mr(dyx, meth="mr_ivw")$pval
}


mr_scatter_plot(mr(dyx), dyx)


hist(rxy_p)
hist(ryx_p)




g1 -> x <- u <- g2
g1 -> y <- u <- g2



