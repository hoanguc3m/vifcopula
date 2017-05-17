library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
t_max = 1000
n_max = 100
datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, family = 1)
datagen <- datagen_gauss
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 50, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_gauss <- vifcopula::vifcop(data,init,other)


pdf(file='img/Gaussian1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_gauss$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_gauss$mean_iv[(t_max+1):(t_max+n_max)] , xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 50, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type
sum(vi_gauss$cop_type == datagen_gauss$family)

################################################################################
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = 2,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = cbind(datagen$family,datagen$family),
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type

################################################################################

datagen_student <- fcopsim(t_max = 1000, n_max = 100, family = 2)
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)

pdf(file='img/Student1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,3))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_student$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_student$mean_iv[seq(from = t_max+1, to = t_max+2*n_max, by = 2 )], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta2, vi_student$mean_iv[seq(from = t_max+2, to = t_max+2*n_max, by = 2 )], xlab = expression(nu[t]), ylab = expression(nu[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, family = 3)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)

pdf(file='img/Clayton1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_clayton$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_clayton$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = gauss_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_clayton <- vifcopula::vifcop(data,init,other)
vi_clayton$cop_type
sum(vi_clayton$cop_type == datagen_clayton$family)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, family = 4)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)

pdf(file='img/Gumbel1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_gumbel$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_gumbel$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = gauss_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gumbel <- vifcopula::vifcop(data,init,other)
vi_gumbel$cop_type
sum(vi_gumbel$cop_type == datagen_gumbel$family)
################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, family = 5)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_frank <- vifcopula::vifcop(data,init,other)

pdf(file='img/Frank1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_frank$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_frank$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = gauss_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_frank <- vifcopula::vifcop(data,init,other)
vi_frank$cop_type
sum(vi_frank$cop_type == datagen_frank$family)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, family = 6)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)

pdf(file='img/Joe1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_joe$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_joe$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = gauss_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_joe <- vifcopula::vifcop(data,init,other)
vi_joe$cop_type
sum(vi_joe$cop_type == datagen_joe$family)

init <- list(copula_type = matrix(1,n_max,1),
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_joe <- vifcopula::vifcop(data,init,other)
vi_joe$cop_type
sum(vi_joe$cop_type == datagen_joe$family)


###############################################################################

copfamily = matrix(sample(c(1,3,4,5,6),size = 100, replace = T),ncol=1)
datagen_mix <- fcopsim(t_max = 1000, n_max = 100, family = copfamily)
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = gauss_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix <- vifcopula::vifcop(data,init,other)

sum(vi_mix$cop_type == datagen_mix$family)

pdf(file='img/Mix1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_mix$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_mix$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()
