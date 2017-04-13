library(TwoSampleMR)
source("simulation-functions.r")


parameters <- init_parameters(nid=10000, nsnp_x=10, var_gx.x=0.1, var_x.y=0.1^2, nsnp_y=10, var_gy.y=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- add_u(parameters, nsnp_u=10, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
parameters <- generate_system_effects(parameters)
sim1 <- simulate_population(parameters)
sim2 <- simulate_population(parameters)
ef1 <- system_effs(sim1)
ef2 <- system_effs(sim2)
d <- system_dat(ef1$x, ef2$y)

cor(sim1$y, sim1$x)
mr(filter(d, inst == "x", pval.exposure < 5e-8))

