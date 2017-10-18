library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
t_max = 1000
n_max = 100
k_max = 6
gid = sample(1:(k_max-1),n_max,replace = T)

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
gauss_latent_init <- matrix(1, nrow = k_max-1, ncol = 1)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max,
                         family = 1, family_latent = 1,
                         gid = gid, structfactor = 3, seed_num = 6)


datagen <- datagen_gauss
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = 4)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)

vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$structfactor <- 3
comparefcop(datagen, vi_gauss)
tail(vi_gauss$mean_iv,105)

plot(datagen$theta, tail(vi_gauss$mean_iv,100), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
text(datagen$theta, tail(vi_gauss$mean_iv,100), label = gid)
abline(a= 0, b=1, col="red")

plot(c(datagen_gauss$theta_latent, datagen_gauss$theta),
     tail(vi_gauss$mean_iv,n_max+k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
text(datagen$theta, tail(vi_gauss$mean_iv,100), label = gid)
abline(a= 0, b=1, col="red")



other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 50, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type
sum(vi_gauss$cop_type == datagen_gauss$family)
sum(vi_gauss$latent_copula_type == datagen_gauss$family_latent)
sum(vi_gauss$gid == datagen_gauss$gid)

###################
indep_init <- matrix(1, nrow = n_max, ncol = 1)
indep_init[10] <- indep_init[20] <- 0
init <- list(copula_type = indep_init,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
vi_gauss <- vifcopula::vifcop(data,init,other)

#################################################################################
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_latent = 3, seed_num = 10,
                           structfactor = 3)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)
plot(datagen_clayton$theta, tail(vi_clayton$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

plot(c(datagen_clayton$theta_latent, datagen_clayton$theta),
     tail(vi_clayton$mean_iv,n_max+k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_clayton <- vifcopula::vifcop(data,init,other)

sum(vi_clayton$cop_type == datagen$family)
sum(vi_clayton$latent_copula_type == datagen$family_latent)


#################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, family_latent = 4, seed_num = 100,
                          structfactor = 3)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)

plot(datagen$theta, tail(vi_gumbel$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

plot(c(datagen_gumbel$theta_latent, datagen_gumbel$theta),
     tail(vi_gumbel$mean_iv,n_max+k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gumbel <- vifcopula::vifcop(data,init,other)

sum(vi_gumbel$cop_type == datagen$family)
sum(vi_gumbel$latent_copula_type == datagen$family_latent)

plot(datagen$theta, tail(vi_gumbel$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
plot(datagen$theta_latent, head(tail(vi_gumbel$mean_iv,n_max+k_max-1),k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

#################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_latent = 5, seed_num = 100,
                         structfactor = 3)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_frank <- vifcopula::vifcop(data,init,other)

plot(datagen$theta, tail(vi_frank$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

plot(c(datagen_frank$theta_latent, datagen_frank$theta),
     tail(vi_frank$mean_iv,n_max+k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")



init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_frank <- vifcopula::vifcop(data,init,other)
sum(vi_frank$cop_type == datagen$family)
sum(vi_frank$latent_copula_type == datagen$family_latent)


#################################################################################
datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 6, family_latent = 6,  seed_num = 100,
                       structfactor = 3)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)

plot(datagen$theta, tail(vi_joe$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


plot(c(datagen_joe$theta_latent, datagen_joe$theta),
     tail(vi_joe$mean_iv,n_max+k_max-1), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_joe <- vifcopula::vifcop(data,init,other)

sum(vi_joe$cop_type == datagen$family)
sum(vi_joe$latent_copula_type == datagen$family_latent)
#################################################################################
copfamily = matrix(sample(c(1,2,3,4,5,6),size = n_max, replace = T),ncol=1)
latentcopfamily = matrix(sample(c(1,2,3,4,5,6),size = k_max-1, replace = T),ncol=1)

copfamily1 = matrix(sample(c(1,2,3,4,5,6),size = n_max, replace = T),ncol=1)
latentcopfamily1 = matrix(sample(c(1,2,3,4,5,6),size = k_max-1, replace = T),ncol=1)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_latent = latentcopfamily,
                       seed_num = 100, structfactor = 3)
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix <- vifcopula::vifcop(data,init,other)


plot(datagen$theta, tail(vi_mix$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

init <- list(copula_type = copfamily1,
             latent_copula_type = latentcopfamily1,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix <- vifcopula::vifcop(data,init,other)

sum(vi_mix$cop_type == datagen$family)
sum(vi_mix$latent_copula_type == datagen$family_latent)

#################################################################################
#################################################################################
#################################################################################
datagen_Student <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 2, family_latent = 2,  seed_num = 100,
                       structfactor = 3)
datagen <- datagen_Student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)

plot(datagen$theta, abs(tail(vi_student$mean_iv,2*n_max)[seq(1,200, by = 2)]), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
plot(datagen$theta2, tail(vi_student$mean_iv,2*n_max)[seq(2,200, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
plot(datagen$theta_latent, vi_student$mean_iv[seq(t_max*k_max+1,t_max*k_max+k_max*2-2,by = 2)])
plot(datagen$theta2_latent, vi_student$mean_iv[seq(t_max*k_max+2,t_max*k_max+k_max*2-1,by = 2)])


abline(a= 0, b=1, col="red")

init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_student <- vifcopula::vifcop(data,init,other)

sum(vi_student$cop_type == datagen$family)
sum(vi_student$latent_copula_type == datagen$family_latent)
#################################################################################
