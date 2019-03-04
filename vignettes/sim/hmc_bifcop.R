setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
set.seed(0)
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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gauss <- vifcopula::hmcfcop(data,init,other)
plot(hmc_gauss)
compare_vi_hmc(vi_gauss, hmc_gauss)
################################################################################

datagen_student <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                           family = 2, seed_num = 100,
                           structfactor = 2, tau_range = c(0.2,0.8),
                           tau_latent_range = c(0.2,0.8), family_latent = sample(c(1,3,4,5,6),size = n_max, replace = T))
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_student <- vifcopula::hmcfcop(data,init,other)

################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_latent = 3, seed_num = 0,
                           tau_range = c(0.2,0.7), tau_latent_range = c(0.2,0.7),
                           structfactor = 2)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_clayton <- vifcopula::hmcfcop(data,init,other)

################################################################################
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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gumbel <- vifcopula::hmcfcop(data,init,other)

################################################################################

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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_frank <- vifcopula::hmcfcop(data,init,other)

###############################################################################

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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_joe <- vifcopula::hmcfcop(data,init,other)

###############################################################################

copfamily = sample(c(1,2,3,4,5,6), size = n_max, replace = T)
latentcopfamily = sample(c(1,3,4,5,6),size = n_max, replace = T)

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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_mix <- vifcopula::hmcfcop(data,init,other)


###############################################################################

time_vi <- c(vi_gauss$time, vi_student$time, vi_clayton$time, vi_gumbel$time, vi_frank$time, vi_joe$time, vi_mix$time)
time_hmc <- c(hmc_gauss$time, hmc_student$time, hmc_clayton$time, hmc_gumbel$time, hmc_frank$time, hmc_joe$time, hmc_mix$time)

time_tab <- rbind(time_vi, time_hmc)
print(xtable(time_tab, digits = 0))

compare_vi_hmc(vi_gauss, hmc_gauss)
compare_vi_hmc(vi_student, hmc_student)
compare_vi_hmc(vi_clayton, hmc_clayton)
compare_vi_hmc(vi_gumbel, hmc_gumbel)
compare_vi_hmc(vi_frank, hmc_frank)
compare_vi_hmc(vi_joe, hmc_joe)
compare_vi_hmc(vi_mix, hmc_mix)

#############################################################################

pdf(file='img/VIHMCbifcop.pdf', width = 18, height = 12)
par(mfrow =c(4,6))
par(mar=c(5,5,3,1))

plot(get_v0_sd(hmc_gauss), get_v0_sd(vi_gauss),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_student), get_v0_sd(vi_student),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_clayton), get_v0_sd(vi_clayton),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_gumbel), get_v0_sd(vi_gumbel),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_frank), get_v0_sd(vi_frank),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_joe), get_v0_sd(vi_joe),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_gauss), get_v_sd(vi_gauss),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_student), get_v_sd(vi_student),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_clayton), get_v_sd(vi_clayton),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_gumbel), get_v_sd(vi_gumbel),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_frank), get_v_sd(vi_frank),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_joe), get_v_sd(vi_joe),
     xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])),
     main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gauss), get_theta_sd(vi_gauss),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_student), get_theta_sd(vi_student),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_clayton), get_theta_sd(vi_clayton),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gumbel), get_theta_sd(vi_gumbel),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_frank), get_theta_sd(vi_frank),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_joe), get_theta_sd(vi_joe),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

plot.new()

plot(get_theta2_sd(hmc_student), get_theta2_sd(vi_student),
     xlab = expression(sd(theta2[hmc]) ), ylab = expression(sd (theta2[vi])),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot.new()
plot.new()
plot.new()
plot.new()
dev.off()

###############################################################################
pdf(file='img/VIvsHMCbifcop.pdf', width = 18, height = 6)
par(mfrow =c(2,6))
par(mar=c(5,5,3,1))

plot(get_theta(hmc_gauss), get_theta(vi_gauss) ,
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta(hmc_student), get_theta(vi_student),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta(hmc_clayton), get_theta(vi_clayton),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta(hmc_gumbel), get_theta(vi_gumbel),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta(hmc_frank), get_theta(vi_frank),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta(hmc_joe), get_theta(vi_joe),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")


plot(get_theta_sd(hmc_gauss), get_theta_sd(vi_gauss),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_student), get_theta_sd(vi_student),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_clayton), get_theta_sd(vi_clayton),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gumbel), get_theta_sd(vi_gumbel),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_frank), get_theta_sd(vi_frank),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_joe), get_theta_sd(vi_joe),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

dev.off()

