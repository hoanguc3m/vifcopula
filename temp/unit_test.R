library("Rcpp")
data <- list(u = u, n_max = n_max, n_group = n_group, t_max = t_max, k_max = k_max, gid = gid, structfactor = 1)
init <- list(copula_type = matrix(copula_type), v = matrix(v0), par = matrix(par))
other <- list(seed = 126, core = 8, iter = 1000, n_monte_carlo_grad = 10, n_monte_carlo_elbo = 100, eval_elbo = 100,
              adapt_bool = F, adapt_val = 1, adapt_iterations = 50, tol_rel_obj = 0.1)

out <- vifcopula::vifcop(data,init,other)
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

cop <- BiCop(family = 5, par = 50)
BiCopDeriv(0.1, 0.5, cop, deriv = "par") /  BiCopPDF(0.1, 0.5, cop)
BiCopDeriv(0.1, 0.5, cop, deriv = "par", log =T)
log(BiCopPDF(0.1, 0.5, cop))
BiCopDeriv(0.1, 0.5, cop, deriv = "u2") /  BiCopPDF(0.1, 0.5, cop)
u1 = runif(1000)
u2 = runif(1000)
sum(BiCopDeriv(u1, u2, cop, deriv = "par") /  BiCopPDF(u1, u2, cop))
sum(log(BiCopPDF(u1, u2, cop)))
#######################################################

datagen <- fcopsim(t_max = 1000, n_max = 100, family = 3)
data <- list(u = datagen$u, n_max = datagen$n_max, n_group = n_group, t_max = datagen$t_max, k_max = datagen$k_max, gid = datagen$gid, structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family, v = datagen$v, par = datagen$theta, par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000, n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10, eval_elbo = 100,
              adapt_bool = F, adapt_val = 1, adapt_iterations = 50, tol_rel_obj = 0.1)
out <- vifcopula::vifcop(data,init,other)

plot(datagen$v, out$mean_iv[1:t_max])
abline(a= 0, b=1, col="red")
plot(datagen$theta, out$mean_iv[(t_max+1):(t_max+n_max)])
abline(a= 0, b=1, col="red")
hist(out$sample_iv[,100])
hist(out$sample_iv[,1050])

BiCopPar2Tau(family = c(3), par = c(50,75,100))
BiCopPar2Tau(family = c(4), par = c(30,50,100))
BiCopPar2Tau(family = c(5), par = c(50,60,100))
BiCopPar2Tau(family = c(6), par = c(25,30,50))
