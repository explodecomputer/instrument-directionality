library(TwoSampleMR)
source("simulation-functions.r")
source("strategies.r")




parameters <- init_parameters(nsnp_x=10, var_gx.x=0.1, var_x.y=0.1^2, nsnp_y=10, var_gy.y=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- generate_system_effects(parameters)
sim1 <- simulate_population(parameters, 10000)
sim2 <- simulate_population(parameters, 10000)
ef1 <- system_effs(sim1)
ef2 <- system_effs(sim2)
d <- system_dat(ef1$x, ef2$y)

cor(sim1$y, sim1$x)
mr(filter(d, inst == "x", pval.exposure < 5e-8))




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


