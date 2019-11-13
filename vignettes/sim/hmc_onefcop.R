setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
t_max = 1000
n_max = 100
gauss_init <- matrix(1, nrow = n_max, ncol = 1)
seed_num = 9876

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, family = 1, seed_num = seed_num)
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
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gauss <- vifcopula::hmcfcop(data,init,other)
vi_gauss <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_student <- fcopsim(t_max = 1000, n_max = 100, family = 2, seed_num = seed_num)
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_student <- vifcopula::hmcfcop(data,init,other)
vi_student <- vifcopula::vifcop(data,init,other)

################################################################################
compare
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, family = 3, seed_num = seed_num)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_clayton <- vifcopula::hmcfcop(data,init,other)
vi_clayton <- vifcopula::vifcop(data,init,other)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, family = 4, seed_num = seed_num)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gumbel <- vifcopula::hmcfcop(data,init,other)
vi_gumbel <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, family = 5, seed_num = seed_num)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_frank <- vifcopula::hmcfcop(data,init,other)
vi_frank <- vifcopula::vifcop(data,init,other)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, family = 6, seed_num = seed_num)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_joe <- vifcopula::hmcfcop(data,init,other)
vi_joe <- vifcopula::vifcop(data,init,other)

###############################################################################

datagen_BB1 <- fcopsim(t_max = 1000, n_max = 100, family = 7, seed_num = seed_num)
datagen <- datagen_BB1
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_BB1 <- vifcopula::hmcfcop(data,init,other)
vi_BB1 <- vifcopula::vifcop(data,init,other)


###############################################################################

copfamily = matrix(sample(c(1,2,3,4,5,6,7),size = 100, replace = T),ncol=1)
datagen_mix <- fcopsim(t_max = 1000, n_max = 100, family = copfamily, family_latent = 0, seed_num = seed_num)
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen_mix$family)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_mix <- vifcopula::hmcfcop(data,init,other)
vi_mix <- vifcopula::vifcop(data,init,other)


###############################################################################
#############################################################################

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

pdf(file='img/VIHMCOnefcop.pdf', width = 19, height = 9)
par(mfrow =c(3,7))
par(mar=c(5,5,3,1))
cex_main = 1.25
cex_lab = 1.5
cex_axis = 0.75

plot(get_v0_sd(hmc_gauss), get_v0_sd(vi_gauss),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gaussian one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_student), get_v0_sd(vi_student),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Student one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_clayton), get_v0_sd(vi_clayton),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Clayton one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_gumbel), get_v0_sd(vi_gumbel),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gumbel one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_frank), get_v0_sd(vi_frank),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Frank one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_joe), get_v0_sd(vi_joe),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Joe one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_BB1), get_v0_sd(vi_BB1),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " BB1 one factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")



plot(get_theta_sd(hmc_gauss), get_theta_sd(vi_gauss),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Gaussian one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_student), get_theta_sd(vi_student),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_clayton), get_theta_sd(vi_clayton),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Clayton one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gumbel), get_theta_sd(vi_gumbel),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Gumbel one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_frank), get_theta_sd(vi_frank),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Frank one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_joe), get_theta_sd(vi_joe),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Joe one factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_BB1), get_theta_sd(vi_BB1),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe one factor copula")
abline(a= 0, b=1, col="red")

plot.new()

plot(get_theta2_sd(hmc_student), get_theta2_sd(vi_student),
     xlab = expression(sd(nu[hmc]) ), ylab = expression(sd (nu[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot.new()
plot.new()
plot.new()
plot.new()
plot(get_theta2_sd(hmc_BB1), get_theta2_sd(vi_BB1),
     xlab = expression(sd(delta[hmc]) ), ylab = expression(sd (delta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
        #main = " BB1 one factor copula")
abline(a= 0, b=1, col="red")

dev.off()
