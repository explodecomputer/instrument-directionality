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

r <- rep(0, 100)
for(i in 1:100)
{
	message(i)
	nid <- 10000
	nsnp <- 20
	g1 <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	g2 <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	x <- makePhen(rep(0.01, nsnp), g1)
	y <- makePhen(c(0.4, rep(0.01, nsnp)), cbind(x, g2))
	d <- get_effs(y, x, g2)
	r[i] <- mr(d, meth="mr_ivw")$pval
}
hist(r)

r1 <- rep(0, 100)
r2 <- rep(0, 100)
for(i in 1:100)
{
	message(i)
	nid <- 10000
	nsnp <- 20
	g1 <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	g2 <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	x <- makePhen(rep(sqrt(0.01), nsnp), g1)
	y <- makePhen(c(sqrt(0.3), rep(sqrt(0.01), nsnp)), cbind(x, g2))
	d <- get_effs(x, y, g1)
	r1[i] <- mr(d, meth="mr_ivw")$pval
	d <- get_effs(y, x, g2)
	r2[i] <- mr(d, meth="mr_ivw")$pval
}

hist(r1)
hist(r2, 20)

sim <- function()
{
	nid <- 100000
	nsnp <- 20
	eff_gx <- rnorm(nsnp)
	eff_gy <- rnorm(nsnp)
	eff_xy <- 0.3
	vgx <- 0.2
	vgy <- 0.2

	Gx <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	Gy <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)
	gx <- Gx %*% eff_gx
	gy <- Gy %*% eff_gy

	x <- scale(gx) * sqrt(vgx) + rnorm(nid, sd=sqrt(1-vgx))
	y <- scale(gy) * sqrt(vgy) + scale(x) * sqrt(eff_xy) + rnorm(nid, sd=sqrt(1-vgx-eff_xy))

	# Get pop1
	get_effs


}




ninst1 <- 50
ninst2 <- 50
ninstu <- 50
nid1 <- 80000
nid2 <- 80000
nidu <- 80000

# Causal effects
var_xy <- 0.3
var_ux <- 0.3
var_uy <- 0

# Genetic effects
var_g1x <- runif(1, 0.05, 0.3)
var_g2y <- runif(1, 0.05, 0.3)
var_guu <- runif(1, 0.05, 0.3)

# Horizontal pleiotropy
# var_g1y <- min(rbeta(1, 1, 100), 0.2)
var_g1y <- 0
# var_g2x <- min(rbeta(1, 1, 100), 0.2)
var_g2x <- 0
# mu_g1y <- min(rbeta(1, 1, 200), 0.05) * sample(c(1,-1), 1)
# mu_g2x <- min(rbeta(1, 1, 300), 0.05) * sample(c(1,-1), 1)
mu_g1y <- 0
mu_g2x <- 0

# mr_method <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
# a <- run_sim(nid1 = nid1, nid2 = nid2, nidu = nidu, ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x, mr_method = mr_method)

effs <- make_effs(ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x)

pop1 <- make_pop(effs, nid1)
pop2 <- make_pop(effs, nid2)
popu <- make_pop(effs, nidu)
ss1 <- get_summary_stats(pop1, pop2, popu)

datux <- make_datg(ss$gw$u, ss$gw$x) %>% filter(inst == "u")
datuy <- make_datg(ss$gw$u, ss$gw$y) %>% filter(inst == "u")

index <- ss$gw$u$inst == "u"
summary(lm(ss$gw$x$b[index] ~ ss$gw$u$b[index]))
summary(lm(ss$gw$y$b[index] ~ ss$gw$u$b[index] + ss$gw$x$b[index]))



ninst1 <- 50
ninst2 <- 50
ninstu <- 50
nid1 <- 80000
nid2 <- 80000
nidu <- 80000

# Causal effects
var_xy <- 0.3
var_ux <- 0.3
var_uy <- 0.3

# Genetic effects
var_g1x <- runif(1, 0.05, 0.3)
var_g2y <- runif(1, 0.05, 0.3)
var_guu <- runif(1, 0.05, 0.3)

# Horizontal pleiotropy
# var_g1y <- min(rbeta(1, 1, 100), 0.2)
var_g1y <- 0
# var_g2x <- min(rbeta(1, 1, 100), 0.2)
var_g2x <- 0
# mu_g1y <- min(rbeta(1, 1, 200), 0.05) * sample(c(1,-1), 1)
# mu_g2x <- min(rbeta(1, 1, 300), 0.05) * sample(c(1,-1), 1)
mu_g1y <- 0
mu_g2x <- 0

# mr_method <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
# a <- run_sim(nid1 = nid1, nid2 = nid2, nidu = nidu, ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x, mr_method = mr_method)

effs <- make_effs(ninst1 = ninst1, ninst2 = ninst2, ninstu = ninstu, var_xy = var_xy, var_ux = var_ux, var_uy = var_uy, var_g1x = var_g1x, var_g2y = var_g2y, var_guu = var_guu, var_g1y = var_g1y, var_g2x = var_g2x, mu_g1y = mu_g1y, mu_g2x = mu_g2x)

pop1 <- make_pop(effs, nid1)
pop2 <- make_pop(effs, nid2)
popu <- make_pop(effs, nidu)
ss2 <- get_summary_stats(pop1, pop2, popu)

datux <- make_datg(ss$gw$u, ss$gw$x) %>% filter(inst == "u")
datuy <- make_datg(ss$gw$u, ss$gw$y) %>% filter(inst == "u")

index <- ss$gw$u$inst == "u"
summary(lm(ss$gw$x$b[index] ~ ss$gw$u$b[index]))
summary(lm(ss$gw$y$b[index] ~ ss$gw$u$b[index] + ss$gw$x$b[index]))
summary(lm(ss$gw$y$b[index] ~ ss$gw$u$b[index]))
summary(lm(ss$gw$y$b[index] ~ ss$gw$x$b[index]))


datux$r <- residuals(lm(beta.exposure ~ beta.outcome, datux))
datux$y <- ss$gw$y$b[index]

ggplot(datux, aes(x=beta.exposure, y=y)) +
geom_point(aes(colour=inst))
ggplot(datux, aes(x=r, y=y)) +
geom_point(aes(colour=inst))




# Not confounder
dat1ux <- make_datg(ss1$gw$u, ss1$gw$x) %>% filter(inst == "u")
dat1uy <- make_datg(ss1$gw$u, ss1$gw$y) %>% filter(inst == "u")


# Counfounder
dat2ux <- make_datg(ss2$gw$u, ss2$gw$x) %>% filter(inst == "u")
dat2uy <- make_datg(ss2$gw$u, ss2$gw$y) %>% filter(inst == "u")


index <- ss1$gw$u$inst == "u"
summary(lm(ss1$gw$x$b[index] ~ ss1$gw$u$b[index]))
summary(lm(ss1$gw$y$b[index] ~ ss1$gw$u$b[index] + ss1$gw$x$b[index]))
summary(lm(ss1$gw$y$b[index] ~ ss1$gw$u$b[index]))
summary(lm(ss1$gw$y$b[index] ~ ss1$gw$x$b[index]))

index <- ss2$gw$u$inst == "u"
summary(lm(ss2$gw$x$b[index] ~ ss2$gw$u$b[index]))
summary(lm(ss2$gw$y$b[index] ~ ss2$gw$u$b[index] + ss2$gw$x$b[index]))
summary(lm(ss2$gw$y$b[index] ~ ss2$gw$u$b[index]))
summary(lm(ss2$gw$y$b[index] ~ ss2$gw$x$b[index]))


summary(lm(ss2$gw$u$b[index] ~ ss2$gw$y$b[index]))
summary(lm(ss2$gw$u$b[index] ~ ss2$gw$y$b[index] + ss2$gw$x$b[index]))

summary(lm(ss1$gw$u$b[index] ~ 0 + ss1$gw$y$b[index]))
summary(lm(ss1$gw$u$b[index] ~ 0 + ss1$gw$x$b[index]))
summary(lm(ss1$gw$u$b[index] ~ ss1$gw$y$b[index] + ss1$gw$x$b[index]))

cit.cp(L = ss1$gw$y$b[index], G = ss1$gw$x$b[index], T = ss1$gw$u$b[index])

cit.cp(L = ss2$gw$y$b[index], G = ss2$gw$x$b[index], T = ss2$gw$u$b[index])




dat2xy <- make_datg(ss2$gw$x, ss2$gw$y) %>% filter(inst == "u")
mr(dat2xy)
var_xy

m <- matrix(c(1, 0.3, 0.3, 0, 1, 0.09, 0, 0, 1), 3)
solve(t(m))

m <- matrix(c(1, 0.3, 0.3, 0, 1, 0, 0, 0, 1), 3)
solve(m)

