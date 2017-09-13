library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)

t_max = 1000
n_max = 100
k_max = 6
gid = sample(1:(k_max-1),n_max,replace = T)

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
gauss_latent_init <- matrix(1, nrow = n_max, ncol = 1)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 1, family_latent = 1, seed_num = 100,
                         structfactor = 2)
datagen <- datagen_gauss
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
vi_gauss <- vifcopula::vifcop(data,init,other)
tail(vi_gauss$mean_iv,105)

plot(datagen$theta[,2], tail(vi_gauss$mean_iv,100), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 50, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type
sum(vi_gauss$cop_type[,1] == datagen_gauss$family[,1])
sum(vi_gauss$latent_copula_type[,1] == datagen_gauss$family_latent[,1])



#################################################################################
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, seed_num = 0,
                           structfactor = 2)
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

sum(vi_clayton$cop_type[,1] == datagen$family[,1])
sum(vi_clayton$latent_copula_type[,1] == datagen$family_latent[,1])


#################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, seed_num = 100,
                          structfactor = 2)
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

plot(datagen$theta[,2], tail(vi_gumbel$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
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

sum(vi_gumbel$cop_type[,1] == datagen$family[,1])
sum(vi_gumbel$latent_copula_type[,1] == datagen$family_latent[,1])


#################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, seed_num = 100,
                         structfactor = 2)
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

plot(datagen$theta[,2], tail(vi_frank$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
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
sum(vi_frank$cop_type[,1] == datagen$family[,1])
sum(vi_frank$latent_copula_type[,1] == datagen$family_latent[,1])


#################################################################################
datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 6, seed_num = 100,
                       structfactor = 2)
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

plot(datagen$theta[,2], tail(vi_joe$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
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

sum(vi_joe$cop_type[,1] == datagen$family[,1])
sum(vi_joe$latent_copula_type[,1] == datagen$family_latent[,1])
#################################################################################
copfamily = matrix(sample(c(1,3,4,5,6),size = n_max, replace = T),ncol=1)
latentcopfamily = matrix(sample(c(1,3,4,5,6),size = k_max-1, replace = T),ncol=1)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_latent = latentcopfamily,
                       seed_num = 100, structfactor = 2)
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


plot(datagen$theta[,2], tail(vi_mix$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
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
vi_mix <- vifcopula::vifcop(data,init,other)

sum(vi_mix$cop_type[,1] == datagen$family[,1])
sum(vi_mix$latent_copula_type[,1] == datagen$family_latent[,1])

#################################################################################
#################################################################################
