# library(devtools)
# install_github("hoanguc3m/vifcopula")
num_rep = 100

setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100

require(parallel) # one of the core R packages
require(doParallel)
# require(multicore); require(doMC) # alternative to parallel/doParallel
# require(Rmpi); require(doMPI) # to use Rmpi as the back-end
library(foreach)
seed_collection <- sample(1:1000000, num_rep, replace=F)

task_onefcop <- function(seed_num, family){
    set.seed(seed_num)

    #gauss_init <- matrix(family, nrow = n_max, ncol = 1)
    copfamily_init <- sample(c(1,2,3,4,5,6,7),size = 100, replace = T)

    datagen <- fcopsim(t_max = t_max, n_max = n_max, family = family, seed_num = seed_num)
    data <- list(u = datagen$u,
                 n_max = datagen$n_max,
                 t_max = datagen$t_max,
                 k_max = datagen$k_max,
                 gid = datagen$gid,
                 structfactor = datagen$structfactor)
    init <- list(copula_type = datagen$family)
    other <- list(seed = seed_num, core = 8, iter = 1000,
                  n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
                  eval_elbo = 100, adapt_bool = F, adapt_val = 1,
                  adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
    vi <- vifcopula::vifcop(data,init,other)


    init <- list(copula_type = copfamily_init)
    other <- list(seed = seed_num, core = 8, iter = 1000,
                  n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
                  eval_elbo = 100, adapt_bool = F, adapt_val = 1,
                  adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T, modelselect = T)
    vi_rng <- vifcopula::vifcop(data,init,other)

    c( correct = sum(vi_rng$cop_type == datagen$family),
       iteration = vi_rng$iteration,
       time_vi = vi$time,
       time_rng = vi_rng$time,
       vi_criteria = vi$criteria,
       vi_rng_criteria = vi_rng$criteria,
       seed_num = seed_num)

}


nCores <- 4
registerDoParallel(nCores)
# registerDoMC(nCores) # alternative to registerDoParallel
# cl <- startMPIcluster(nCores); registerDoMPI(cl) # when using Rmpi as the

library(doRNG, quietly = TRUE)

Data_Gauss <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 1)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Student <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 2)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Clayton <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 3)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Gumbel <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 4)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Frank <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 5)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Joe <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 6)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_BB1 <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = 7)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}

Data_Mix <- foreach(i = 1:num_rep, .combine= 'cbind', .options.RNG = list(seed = 0)) %dopar% {
    copfamily = sample(c(1,2,3,4,5,6,7),size = 100, replace = T)
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- task_onefcop(seed_collection[i], family = copfamily)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
}




#############################################################################
# correct        iteration          time_vi         time_rng     vi_criteria1     vi_criteria2     vi_criteria3     vi_criteria4
# vi_rng_criteria1  vi_rng_criteria2    vi_rng_criteria3    vi_rng_criteria4
vi_gauss <- apply(Data_Gauss, MARGIN = 1, FUN = mean)
vi_student <- apply(Data_Student, MARGIN = 1, FUN = mean)
vi_clayton <- apply(Data_Clayton, MARGIN = 1, FUN = mean)
vi_gumbel <- apply(Data_Gumbel, MARGIN = 1, FUN = mean)
vi_frank <- apply(Data_Frank, MARGIN = 1, FUN = mean)
vi_joe <- apply(Data_Joe, MARGIN = 1, FUN = mean)
vi_BB1 <- apply(Data_BB1, MARGIN = 1, FUN = mean)
vi_mix <- apply(Data_Mix, MARGIN = 1, FUN = mean)


time <- c(vi_gauss[3], vi_student[3], vi_clayton[3], vi_gumbel[3], vi_frank[3], vi_joe[3], vi_BB1[3], vi_mix[3])
print(time, digits = 0)

ELBO_init <- c(vi_gauss[5], vi_student[5], vi_clayton[5], vi_gumbel[5], vi_frank[5], vi_joe[5], vi_BB1[5], vi_mix[5])
print(ELBO_init, digits = 1)

AIC_init <- c(vi_gauss[6], vi_student[6], vi_clayton[6], vi_gumbel[6], vi_frank[6], vi_joe[6], vi_BB1[6], vi_mix[6])
BIC_init <- c(vi_gauss[7], vi_student[7], vi_clayton[7], vi_gumbel[7], vi_frank[7], vi_joe[7], vi_BB1[7], vi_mix[7])
logP_init <- c(vi_gauss[8], vi_student[8], vi_clayton[8], vi_gumbel[8], vi_frank[8], vi_joe[8], vi_BB1[8], vi_mix[8])

init_tab <- rbind(ELBO_init, AIC_init, BIC_init, logP_init)/t_max
library(xtable)
print(xtable(init_tab, digits = 1))

iter_num <- c(vi_gauss[2], vi_student[2], vi_clayton[2], vi_gumbel[2], vi_frank[2], vi_joe[2], vi_BB1[2], vi_mix[2])
print(iter_num, digits = 0)
iter_num <- round(iter_num)

correct_percent <- c(vi_gauss[1], vi_student[1], vi_clayton[1], vi_gumbel[1], vi_frank[1], vi_joe[1], vi_BB1[1], vi_mix[1])
print(correct_percent, digits = 0)
correct_percent <- round(correct_percent)

time_rng <- c(vi_gauss[4], vi_student[4], vi_clayton[4], vi_gumbel[4], vi_frank[4], vi_joe[4], vi_BB1[4], vi_mix[4])
print(time_rng, digits = 0)

ELBO_rng <- c(vi_gauss[9], vi_student[9], vi_clayton[9], vi_gumbel[9], vi_frank[9], vi_joe[9], vi_BB1[9], vi_mix[9])
AIC_rng <- c(vi_gauss[10], vi_student[10], vi_clayton[10], vi_gumbel[10], vi_frank[10], vi_joe[10], vi_BB1[10], vi_mix[10])
BIC_rng <- c(vi_gauss[11], vi_student[11], vi_clayton[11], vi_gumbel[11], vi_frank[11], vi_joe[11], vi_BB1[11], vi_mix[11])
logP_rng <- c(vi_gauss[12], vi_student[12], vi_clayton[12], vi_gumbel[12], vi_frank[12], vi_joe[12], vi_BB1[12], vi_mix[12])


rng_tab <- rbind(iter_num, correct_percent, rbind(ELBO_rng, AIC_rng, BIC_rng, logP_rng)/t_max)

#print(ELBO_rng, digits = 1)
print(xtable(rng_tab, digits = 1))
print(xtable(rng_tab, digits = 0))



plot(Data_Gauss[8,], Data_Gauss[12,])
abline(a= 0, b= 1)
plot(Data_Student[8,], Data_Student[12,])
abline(a= 0, b= 1)
plot(Data_Clayton[8,], Data_Clayton[12,])
abline(a= 0, b= 1)
plot(Data_Gumbel[8,], Data_Gumbel[12,])
abline(a= 0, b= 1)
plot(Data_Frank[8,], Data_Frank[12,])
abline(a= 0, b= 1)
plot(Data_Joe[8,], Data_Joe[12,])
abline(a= 0, b= 1)
plot(Data_BB1[8,], Data_BB1[12,])
abline(a= 0, b= 1)
plot(Data_Mix[8,], Data_Mix[12,])
abline(a= 0, b= 1)

