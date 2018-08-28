# 1-factor
cpar.frk=c(12.2,3.45,4.47,4.47,5.82); d=5
n=10000
set.seed(123)
frkdat=sim1fact(n,cpar.frk,qcondfrk,"frk")
cat("\nFrank 1-factor MLE: standalone R and then f90\n")
out.frk=ml1fact(nq=25,cpar.frk,frkdat,dfrk,LB=-30,UB=30,prlevel=1,mxiter=100)



data <- list(u = frkdat,
             n_max = d,
             t_max = n,
             k_max = 1,
             gid = c(rep(1, d) ),
             structfactor = 1)
init <- list(copula_type = rep(5,d) )
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)
theta = get_theta(vi_frank)


# 2-factor
d=7; be1=c(.7,.6,.7,.6,.7,.6,.7); be2=c(.4,.4,.4,.4,.3,.3,.3)
cpar1.frk=frk.b2cpar(be1); cpar2.frk=frk.b2cpar(be2)
n=1000
set.seed(123)
frkdat=sim2fact(n,cpar1.frk,cpar2.frk,qcondfrk,qcondfrk,"frk","frk")
cat("\nFrank 2-factor MLE: nlm and then pdhessmin\n")
out.frk=ml2fact(nq=25,c(cpar1.frk,cpar2.frk),frkdat,copname="frank",
                LB=-30,UB=30,prlevel=1,mxiter=100)



data <- list(u = frkdat,
             n_max = d,
             t_max = n,
             k_max = 2,
             gid = c(rep(1, d) ),
             structfactor = 2)
init <- list(copula_type = rep(5,d),
                latent_copula_type = rep(5,d))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)
vi_frank$criteria
theta = c(get_theta(vi_frank), get_latent_theta(vi_frank))


##################################################################

gl=gausslegendre(25)
grsize=c(4,4,3)
d=sum(grsize)
n=1000
set.seed(123)

parbi=c(rep(4,11),rep(6,4),rep(6.5,4),rep(7,3))
udatbi=simbifact(n,grsize,cop=5,parbi)
npar=2*d
dstrfrk=list(data=udatbi,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=0)
nllk=f90str2nllk(parbi,dstrfrk)
        out.frk=f90str2nllk(param,dstrfrk)
outb=pdhessminb(c(rep(2,d),rep(3,d)),f90str2nllk,ifixed=rep(FALSE,npar),dstrfrk,
                LB=rep(0,npar), UB=rep(20,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
out.frk=f90str2nllk(outb$parmin,dstrfrk)


data <- list(u = udatbi,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]) ),
             structfactor = 2)

init <- list(copula_type = rep(5,d),
             latent_copula_type = rep(5,d))
other <- list(seed = 12611, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = T, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)
vi_frank$criteria
param = c(get_theta(vi_frank), get_latent_theta(vi_frank))

################################################################

gl=gausslegendre(25)
grsize=c(4,4,3,4,4,3)
d=sum(grsize)
n=1000
# nested-factor copula
mgrp=length(grsize)
set.seed(123)
parne=c(rep(4,3),rep(6,4),rep(6.5,4),rep(4,3),rep(6,4),rep(6.5,4),rep(7,3),rep(7,3))
udatne=simnestfact(n,grsize,cop=5,parne)
dstrfrk=list(data=udatne,copname="frank",quad=gl,repar=0,grsize=grsize)
npar=mgrp+d
nllk=f90str1nllk(parne,dstrfrk)
nllk$fnval
outn= pdhessminb(rep(3,npar),f90str1nllk, ifixed=rep(FALSE,npar), dstrfrk,
                 LB=rep(0,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
outn$fnval
#out.frk=f90str1nllk(param,dstrfrk)
out.frk$fnval

data <- list(u = udatne,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]),
                     rep(4, grsize[4]),rep(5, grsize[5]),rep(6, grsize[6]) ),
             structfactor = 3)
init <- list(copula_type = rep(5,d),
             latent_copula_type = rep(5,length(grsize)))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)
vi_frank$criteria

param = c(get_latent_theta(vi_frank),get_theta(vi_frank))





datagen_frank <- fcopsim(t_max = 1000, n_max = 22, k_max = 7, gid = gid,
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
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
vi_frank <- vifcopula::vifcop(data,init,other)

