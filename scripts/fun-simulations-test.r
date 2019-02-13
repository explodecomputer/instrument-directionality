library(TwoSampleMR)
library(simulateGP)

ss <- simulateGP::create_system(
    nidx=sample(20000:500000, 1),
    nidy=sample(20000:500000, 1),
    nidu=0,
    nu=sample(0:10, 1),
    na=0,
    nb=0,
    var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
    nsnp_x=sample(1:200, 1),
    nsnp_y=sample(1:200, 1),
    var_gx.x=runif(1, 0.01, 0.1),
    var_gy.y=runif(1, 0.01, 0.1),
    var_gx.y=runif(1, 0.001, 0.01),
    mu_gx.y=runif(1, -0.005, 0.005),
    prop_gx.y=runif(1, 0, 1),
    var_gy.x=runif(1, 0.001, 0.01),
    mu_gy.x=runif(1, -0.005, 0.005),
    prop_gy.x=runif(1, 0, 1)
)

dx <- make_dat(ss$x$x, ss$y$y)
dy <- make_dat(ss$y$y, ss$x$x)

mr_wrapper(subset(dx, pval.exposure < 1e-8))
mr_wrapper(subset(dy, pval.exposure < 1e-8))

# OR


parameters <- init_parameters(
    var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
    nsnp_x=sample(1:200, 1),
    nsnp_y=sample(1:200, 1),
    var_gx.x=runif(1, 0.01, 0.1),
    var_gy.y=runif(1, 0.01, 0.1),
    var_gx.y=runif(1, 0.001, 0.01),
    mu_gx.y=runif(1, -0.005, 0.005),
    prop_gx.y=runif(1, 0, 1),
    var_gy.x=runif(1, 0.001, 0.01),
    mu_gy.x=runif(1, -0.005, 0.005),
    prop_gy.x=runif(1, 0, 1)
)

parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)

parameters <- generate_system_effects(parameters)
sim1 <- simulate_population(parameters, 10000)
sim2 <- simulate_population(parameters, 20000)
ef1 <- system_effs(sim1)
ef2 <- system_effs(sim2)
d <- make_dat(ef1$x, ef2$y)


mr_wrapper(subset(dx, pval.exposure < 1e-8))
mr_wrapper(subset(dy, pval.exposure < 1e-8))
