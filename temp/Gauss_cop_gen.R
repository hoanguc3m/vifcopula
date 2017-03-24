set.seed(126)
library("VineCopula")
#library("statmod")


n_max  <- 10           # Set the number of time series
n_group <- 10           # Set the number of G groups
t_max  <-  100         # Set the number of periods
k_max  <-  1         # Set the number of k latents

gid <- sample(1:n_group, n_max, replace = T)
copula_type <- rep(1,n_max)
rho0  <-  runif(n_max,min = 0.1, max = 0.9)
burn_in  <-  0
nsim <- burn_in + t_max

v0  <- runif(nsim)                           # The latent variable
u <- matrix(0,nrow=t_max, ncol = n_max)

for (i in 1:n_max){
    obj <- BiCop(family = 1, par = rho0[i])
    u[,i] <- BiCopCondSim(t_max, cond.val = v0, cond.var = 1, obj)
}
v0_true <- v0
rho0_true <- rho0

par  <-  runif(n_max,min = 0.1, max = 0.9)
v0  <- runif(nsim)                           # The latent variable

data <- list(u = matrix(u,nrow = t_max), n_max = n_max, n_group = n_group, t_max = t_max, k_max = k_max, gid = gid, structfactor = 1)
init <- list(copula_type = matrix(copula_type), v = matrix(v0), par = matrix(par))
other <- list(seed = 126, core = 8, iter = 1000, n_monte_carlo_grad = 5, n_monte_carlo_elbo = 5)

################################################################
# Gen Clayton copula
################################################################

datagen <- fcopsim(t_max = 1000, n_max = 100, family = 3)


