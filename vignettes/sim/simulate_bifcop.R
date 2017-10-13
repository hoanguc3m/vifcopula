library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)

t_max = 1000
n_max = 100
k_max = 6
gid = sample(1:(k_max-1),n_max,replace = T)
copfamily_rng = sample(c(1,3,4,5,6), size = n_max, replace = T)
latentcopfamily_rng = sample(c(1,3,4,5,6),size = n_max, replace = T)

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
comparefcop(datagen, vi_gauss)
#tail(vi_gauss$mean_iv,105)

pdf(file='img/GaussianBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_gauss$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_gauss$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

hist(get_v0(vi_gauss))
hist(get_v(vi_gauss))

plot(datagen$theta_latent, get_latent_theta(vi_gauss), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_gauss), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


init <- list(copula_type = copfamily_rng,
             latent_copula_type = latentcopfamily_rng,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)


other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type
sum(vi_gauss$cop_type == datagen_gauss$family)
sum(vi_gauss$latent_copula_type == datagen_gauss$family_latent)

###################
indep_init <- matrix(1, nrow = n_max, ncol = 1)
indep_init[10] <- indep_init[20] <- 0
init <- list(copula_type = indep_init,
             latent_copula_type = datagen$family_latent,
             v = datagen$v,
             par = datagen$theta,
             par2 = datagen$theta2)
vi_gauss <- vifcopula::vifcop(data,init,other)
num_param(vi_gauss)
#################################################################################


datagen_student <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 2, family_latent = 2, seed_num = 100,
                         structfactor = 2)
datagen <- datagen_student
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


pdf(file='img/StudentBifc.pdf', width = 9, height = 9)
par(mfrow =c(2,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_student$mean_iv,n_max)[seq(1,2*n_max, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta2_latent, tail(vi_student$mean_iv,n_max)[seq(2,2*n_max, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")


plot(datagen$theta, vi_student$mean_iv[seq(t_max*k_max+1,t_max*k_max+2*n_max, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta2, vi_student$mean_iv[seq(t_max*k_max+2,t_max*k_max+2*n_max, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")

dev.off()






other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
vi_gauss$cop_type
sum(vi_gauss$cop_type == datagen_gauss$family)
sum(vi_gauss$latent_copula_type == datagen_gauss$family_latent)



#################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_latent = 3, seed_num = 0,
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

pdf(file='img/ClaytonBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_clayton$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_clayton$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()


init <- list(copula_type = copfamily_rng,
            latent_copula_type = latentcopfamily_rng,
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

pdf(file='img/GumbelBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_gumbel$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_gumbel$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()


init <- list(copula_type = copfamily_rng,
             latent_copula_type = latentcopfamily_rng,
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


#################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_latent = 5, seed_num = 100,
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

pdf(file='img/FrankBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_frank$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_frank$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = copfamily_rng,
             latent_copula_type = latentcopfamily_rng,
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
                       family = 6, family_latent = 6, seed_num = 100,
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

pdf(file='img/JoeBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_joe$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_joe$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

init <- list(copula_type = copfamily,
             latent_copula_type = latentcopfamily,
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
copfamily = sample(c(1,3,4,5,6), size = n_max, replace = T)
latentcopfamily = sample(c(1,3,4,5,6),size = n_max, replace = T)

copfamily1 = sample(c(1,3,4,5,6), size = n_max, replace = T)
latentcopfamily1 = sample(c(0,1,3,4,5,6),size = n_max, replace = T)

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
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)

pdf(file='img/MixBifc.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$theta_latent, tail(vi_mix$mean_iv,n_max), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_mix$mean_iv[(t_max*k_max+1):(t_max*k_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()

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
