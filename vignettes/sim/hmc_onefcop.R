library(vifcopula)
t_max = 1000
n_max = 100
gauss_init <- matrix(1, nrow = n_max, ncol = 1)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, family = 1, seed_num = 0)
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
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_gauss <- vifcopula::hmcfcop(data,init,other)
################################################################################

datagen_student <- fcopsim(t_max = 1000, n_max = 100, family = 2, seed_num = 0)
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_student <- vifcopula::hmcfcop(data,init,other)

################################################################################
compare
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, family = 3, seed_num = 0)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_clayton <- vifcopula::hmcfcop(data,init,other)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, family = 4, seed_num = 0)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_gumbel <- vifcopula::hmcfcop(data,init,other)

################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, family = 5, seed_num = 0)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_frank <- vifcopula::hmcfcop(data,init,other)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, family = 6, seed_num = 0)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_joe <- vifcopula::hmcfcop(data,init,other)

###############################################################################

copfamily = matrix(sample(c(1,2,3,4,5,6),size = 100, replace = T),ncol=1)
datagen_mix <- fcopsim(t_max = 1000, n_max = 100, family = copfamily, family_latent = 0, seed_num = 0)
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen_mix$family)
other <- list(seed = 126, core = 8, iter = 1000, copselect = F)
hmc_mix <- vifcopula::hmcfcop(data,init,other)


###############################################################################
