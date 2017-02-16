
data <- list(u = u, n_max = n_max, n_group = n_group, t_max = t_max, k_max = k_max, gid = gid, structfactor = 1)
init <- list(copula_type = matrix(copula_type), v = matrix(v0), par = matrix(par))
other <- list(seed = 126, core = 8, iter = 1000, n_monte_carlo_grad = 5, n_monte_carlo_elbo = 5)

vifcopula::vifcop(data,init,other)
vifcop(data,init,other)
sourceCpp("src/vifcop.cpp")
sourceCpp("src/logBifcop.cpp")
sourceCpp("src/bicop_independence_log.cpp")
sourceCpp("src/distribution/bicop_normal_log.cpp")

cop <- BiCop(family = 1, par = .5)
cop <- BiCop(family = 2, par = .5, par2 = 5)
sum(log(BiCopPDF(u[,1], v0, cop)))


data <- list(u = u[1,1], n_max = 1, n_group = 1, t_max = 1, k_max = 1, gid = 1, structfactor = 1)
init <- list(copula_type = matrix(rep(2,n_max)), v = matrix(v0), par = matrix(par))
