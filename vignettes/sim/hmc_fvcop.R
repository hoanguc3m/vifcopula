library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100
k_max = 1
g_max = 10

gid = sample(1:g_max,n_max,replace = T)

copfamily_rng = sample(c(1,3,4,5,6), size = n_max, replace = T)
vinecopfamily_rng = sample(c(1,3,4,5,6),size = n_max, replace = T)

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
gauss_vine_init <- matrix(1, nrow = n_max - g_max, ncol = 1)

################################################################################
datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 1, family_vine = 1, seed_num = 100,
                         structfactor = 4)
datagen <- datagen_gauss

data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gauss <- vifcopula::hmcfcop(data,init,other)
plot(hmc_gauss)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_gauss <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_student <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                           family = 2, seed_num = 100,
                           structfactor = 4,
                           family_vine = sample(c(1,3,4,5,6),size = n_max - g_max, replace = T))
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_student <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_student <- vifcopula::vifcop(data,init,other)

################################################################################
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           tau_range = c(0.2,0.8), tau_latent_range = c(0.2,0.5),
                           family = 3, family_vine = 3, seed_num = 84210,
                           structfactor = 4)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_clayton <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, family_vine = 4, seed_num = 100,
                          structfactor = 4)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gumbel <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_vine = 5, seed_num = 100,
                         structfactor = 4)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_frank <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_frank <- vifcopula::vifcop(data,init,other)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       tau_range = c(0.2,0.8), tau_latent_range = c(0.2,0.5),
                       family = 6, family_vine = 6, seed_num = 985412,
                       structfactor = 4)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_joe <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)

###############################################################################

copfamily = sample(c(1,2,3,4,5,6), size = n_max, replace = T)
latentcopfamily = sample(c(1,3,4,5,6),size = n_max - max(gid), replace = T)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_vine = latentcopfamily,
                       seed_num = 100, structfactor = 4)

datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_mix <- vifcopula::hmcfcop(data,init,other)
other <- list(seed = 126, core = 8, iter = 1000,
              tol_rel_obj = 0.01, copselect = F, modelselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)


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
setwd("/home/hoanguc3m/Dropbox/Wp4/")
pdf(file='img/VIHMCfvcop.pdf', width = 21, height = 12)
par(mfrow =c(4,7))
par(mar=c(5,5,3,1))

plot(get_v0(hmc_gauss), get_v0(vi_gauss),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Gaussian truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_student), get_v0(vi_student),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Student truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_clayton), get_v0(vi_clayton),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Clayton truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_gumbel), get_v0(vi_gumbel),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Gumbel truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_frank), get_v0(vi_frank),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Frank truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_joe), get_v0(vi_joe),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Joe truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(get_v0(hmc_mix), get_v0(vi_mix),
     xlab = expression(v0[hmc] ), ylab = expression(v0[vi]),
     main = " Mix truncated factor vine copula")
abline(a= 0, b=1, col="red")


##############################


plot(get_v0_sd(hmc_gauss), get_v0_sd(vi_gauss),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_student), get_v0_sd(vi_student),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_clayton), get_v0_sd(vi_clayton),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_gumbel), get_v0_sd(vi_gumbel),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_frank), get_v0_sd(vi_frank),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_joe), get_v0_sd(vi_joe),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_mix), get_v0_sd(vi_mix),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd(v0[vi])),
     main = "")
abline(a= 0, b=1, col="red")
#########################################


plot( c(get_theta(hmc_gauss),
        get_vine_theta(hmc_gauss)),
      c(get_theta(vi_gauss),
        get_vine_theta(vi_gauss)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_student),
        get_theta2(hmc_student),
        get_vine_theta(hmc_student)),
      c(get_theta(vi_student),
        get_theta2(vi_student),
        get_vine_theta(vi_student)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_clayton),
        get_vine_theta(hmc_clayton)),
      c(get_theta(vi_clayton),
        get_vine_theta(vi_clayton)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_gumbel),
        get_vine_theta(hmc_gumbel)),
      c(get_theta(vi_gumbel),
        get_vine_theta(vi_gumbel)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_frank),
        get_vine_theta(hmc_frank)),
      c(get_theta(vi_frank),
        get_vine_theta(vi_frank)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_joe),
        get_vine_theta(hmc_joe)),
      c(get_theta(vi_joe),
        get_vine_theta(vi_joe)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta(hmc_mix),
        get_vine_theta(hmc_mix)),
      c(get_theta(vi_mix),
        get_vine_theta(vi_mix)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

##################################
plot( c(get_theta_sd(hmc_gauss),
        get_vine_theta_sd(hmc_gauss)),
      c(get_theta_sd(vi_gauss),
        get_vine_theta_sd(vi_gauss)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_student),
        get_theta2_sd(hmc_student),
        get_vine_theta_sd(hmc_student)),
      c(get_theta_sd(vi_student),
        get_theta2_sd(vi_student),
        get_vine_theta_sd(vi_student)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_clayton),
        get_vine_theta_sd(hmc_clayton)),
      c(get_theta_sd(vi_clayton),
        get_vine_theta_sd(vi_clayton)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_gumbel),
        get_vine_theta_sd(hmc_gumbel)),
      c(get_theta_sd(vi_gumbel),
        get_vine_theta_sd(vi_gumbel)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_frank),
        get_vine_theta_sd(hmc_frank)),
      c(get_theta_sd(vi_frank),
        get_vine_theta_sd(vi_frank)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_joe),
        get_vine_theta_sd(hmc_joe)),
      c(get_theta_sd(vi_joe),
        get_vine_theta_sd(vi_joe)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

plot( c(get_theta_sd(hmc_mix),
        get_vine_theta_sd(hmc_mix)),
      c(get_theta_sd(vi_mix),
        get_vine_theta_sd(vi_mix)),
      xlab = expression(theta[hmc] ), ylab = expression(theta[vi]),
      main = "")
abline(a= 0, b=1, col="red")

##################################
dev.off()

