library("Rcpp")
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
cop <- BiCop(family = 3, par = .5)
cop <- BiCop(family = 4, par = 2)
cop <- BiCop(family = 5, par = 2)
cop <- BiCop(family = 6, par = 3)
sum(log(BiCopPDF(u[,1], v0, cop)))


data <- list(u = u[1,1], n_max = 1, n_group = 1, t_max = 1, k_max = 1, gid = 1, structfactor = 1)
init <- list(copula_type = matrix(rep(6,n_max)), v = matrix(v0), par = matrix(par))

r <- stanc(file = "/home/hoanguc3m/NetBeansProjects/stan-develop/src/test/test-models/good/variational/hier_logistic_cp.stan", model_name = "hier_logistic_cp_model")
str(r)
r$model_code
write(r$cppcode, file = "11212121212.txt")


library("VineCopula")

# cop <- BiCop(family = 1, par = 0.5)
# log(BiCopPDF(0.1, 0.1, cop))
#
# BiCopDeriv(u1, u2, cop, deriv = "u2", log = TRUE)
# BiCopDeriv(0.1, 0.5, cop, deriv = "par", log = TRUE)
#
# log(BiCopPDF(0.1, 0.5, cop))
# BiCopDeriv(0.1, 0.1, cop, deriv = "u2") / BiCopPDF(0.1, 0.1, cop)

cop <- BiCop(family = 2, par = 0.5, par2 = 4)
cop <- BiCop(family = 2, par = 0.5, par2 = 6)
cop <- BiCop(family = 2, par = 0.6, par2 = 5)

log(BiCopPDF(0.1, 0.1, cop))

BiCopDeriv(0.1, 0.1, cop, deriv = "u2") / BiCopPDF(0.1, 0.1, cop)
BiCopDeriv(0.1, 0.5, cop, deriv = "par") / BiCopPDF(0.1, 0.5, cop)
BiCopDeriv(0.1, 0.8, cop, deriv = "par2") / BiCopPDF(0.1, 0.8, cop)

cop1 <- BiCop(family = 2, par = 0.6, par2 = 5+0.0001)
cop2 <- BiCop(family = 2, par = 0.6, par2 = 5-0.0001)
(log(BiCopPDF(0.1, 0.9, cop1) ) - log(BiCopPDF(0.1, 0.9, cop2) ))/(0.0002)

sourceCpp("inst/include/unitTests/temp.cpp")

cop <- BiCop(family = 2, par = 0.8, par2 = 5)
BiCopDeriv(0.6, 0.4, cop, deriv = "par2", log = T)
log(BiCopPDF(0.6, 0.4, cop))
BiCopDeriv(0.6, 0.4, cop, deriv = "par2") / BiCopPDF(0.6, 0.4, cop)

cop1 <- BiCop(family = 2, par = 0.8, par2 = 5+0.0001)
cop2 <- BiCop(family = 2, par = 0.8, par2 = 5-0.0001)
(log(BiCopPDF(0.6, 0.4, cop1) ) - log(BiCopPDF(0.6, 0.4, cop2) ))/(0.0002)

cop <- BiCop(family = 4, par = 2)
BiCopDeriv(0.1, 0.1, cop, deriv = "par", log = T)
log(BiCopPDF(0.1, 0.2, cop))
BiCopDeriv(0.1, 0.1, cop, deriv = "par") / BiCopPDF(0.1, 0.1, cop)

BiCopDeriv(0.1, 0.1, cop, deriv = "u1") /  BiCopPDF(0.1, 0.1, cop)

log(BiCopPDF(0.1, 0.3, cop))
BiCopDeriv(0.1, 0.3, cop, deriv = "u2") /  BiCopPDF(0.1, 0.3, cop)

