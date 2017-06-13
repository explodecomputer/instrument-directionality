library(TwoSampleMR)
source("simulation-functions.r")
source("strategies.r")




parameters <- init_parameters(nsnp_x=100, var_gx.x=0.3, var_x.y=0.3^2, nsnp_y=10, var_gy.y=0.1, mu_gx.y = -0.025, prop_gx.y=0.5)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- generate_system_effects(parameters)
sim1 <- simulate_population(parameters, 100000)
sim2 <- simulate_population(parameters, 100000)
ef1 <- system_effs(sim1)
ef2 <- system_effs(sim2)
d <- system_dat(ef1$x, ef2$y)
d1 <- filter(d, inst == "x", pval.exposure < 5e-8)
cor(sim1$y, sim1$x)
m <- mr(d1)

mr_scatter_plot(m, d1)[[1]] + ggplot2::xlim(c(0, max(d1$beta.exposure)))




out <- create_system(nidx=50000, nidy=50000, nidu=100, nu=5, na=5, nb=5, nsnp_x=100, var_gx.x=0.3, var_x.y=0.3^2, nsnp_y=10, var_gy.y=0.1, mu_gx.y = -0.025, prop_gx.y=0.5, simulate_nonxy=FALSE)


# What is the output
# out$x$x is the genetic effects of all SNPs on x trait, using the discovery population for x
# out$y$x is the genetic effects of all SNPs on x trait, using the discovery population for y
# out$y$u1 is the genetic effects of all SNPs on u1 trait, using the discovery population for y


######

ss <- create_system(40000, 50000, 2000, 5, 5, 5, 0.1)
p <- ss$parameters
names(p$u[[1]])

p$u[[15]]

sapply(ss$parameters$u, function(x) x$name_u)

ss1 <- strategies_oracle(ss)
ss2 <- strategies_tophits(ss)
ss3 <- strategies_directionality(ss)


res1xy <- mr_all(ss1$xy)
res1yx <- mr_all(ss1$yx)
res2xy <- mr_all(ss2$xy)
res2yx <- mr_all(ss2$yx)
res3xy <- mr_all(ss3$xy)
res3yx <- mr_all(ss3$yx)


mr_rucker(a$xy)
mr_rucker(a$yx)

mr_median(a$xy)
mr_mode(a$xy)


plot(a$xy$beta.outcome ~ a$xy$beta.exposure)

ss1 <- create_system(400, 500, 200, 5, 5, 5, 0.1)


