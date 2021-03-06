#' @export
compare_sim_vi <- function(datagen,vi){

    if (! all.equal(vi$structfactor, datagen$structfactor))
        stop("Factor copula model needs to be the same ")

    op <- par('bg', 'col', 'col.axis')
    on.exit(par(op))
    v0_vi = get_v0(vi)
    v_vi = get_v(vi)


    theta_vi = get_theta(vi)
    theta2_vi = get_theta2(vi)

    tau_vi = BiCopPar2Tau(family = vi$cop_type, par = theta_vi, par2 = theta2_vi)
    tau = BiCopPar2Tau(family = datagen$family, par = datagen$theta, par2 = datagen$theta2)


    if (vi$structfactor == 1) {
        par(mfrow =c(1,3))
    } else {
        par(mfrow =c(2,3))
    }

    par(mar=c(5,5,3,1))

    plot(datagen$v[,1], v0_vi, xlab = expression(v[gen]), ylab = expression(v[vi]))
    abline(a= 0, b=1, col="red")

    plot(tau, tau_vi, xlab = expression(tau[gen]), ylab = expression(tau[vi]))
    abline(a= 0, b=1, col="red")

    plot(datagen$theta, theta_vi, xlab = expression(theta[gen]), ylab = expression(theta[vi]))
    abline(a= 0, b=1, col="red")

    if (vi$structfactor == 2 | vi$structfactor == 3){
        latent_theta_vi = get_latent_theta(vi)
        latent_theta2_vi = get_latent_theta2(vi)

        plot(datagen$v[,2:vi$k_max], v_vi, xlab = expression(v[gen]), ylab = expression(v[vi]))
        abline(a= 0, b=1, col="red")
        latent_tau_vi = BiCopPar2Tau(family = vi$latent_copula_type,
                                     par = latent_theta_vi, par2 = latent_theta2_vi)
        latent_tau = BiCopPar2Tau(family = datagen$family_latent,
                                  par = datagen$theta_latent, par2 = datagen$theta2_latent)
        plot(latent_tau, latent_tau_vi, xlab = expression(tau_latent[gen]), ylab = expression(tau_latent[vi]))
        abline(a= 0, b=1, col="red")

        plot(datagen$theta_latent, latent_theta_vi, xlab = expression(theta_latent[gen]), ylab = expression(theta_latent[vi]))
        abline(a= 0, b=1, col="red")
    }

    if (vi$structfactor == 4) {
        plot(datagen$theta2, theta2_vi, xlab = expression(theta2[gen]), ylab = expression(theta2[vi]))
        abline(a= 0, b=1, col="red")

        vine_theta_vi = get_vine_theta(vi)
        vine_theta2_vi = get_vine_theta2(vi)

        vine_tau_vi = BiCopPar2Tau(family = vi$vine_copula_type,
                                   par = vine_theta_vi, par2 = vine_theta2_vi)
        vine_tau = BiCopPar2Tau(family = datagen$family_vine,
                                par = datagen$theta_vine, par2 = datagen$theta2_vine)
        plot(vine_tau, vine_tau_vi, xlab = expression(tau_vine[gen]), ylab = expression(tau_vine[vi]))
        abline(a= 0, b=1, col="red")

        plot(datagen$theta_vine, vine_theta_vi, xlab = expression(theta_vine[gen]), ylab = expression(theta_vine[vi]))
        abline(a= 0, b=1, col="red")

        }
    par(mfrow =c(1,1))
}

#' @export
plot.vifcop <- function(vi) {
    op <- par('bg', 'col', 'col.axis')
    on.exit(par(op))
    v0_vi = get_v0(vi)
    v_vi = get_v(vi)


    theta_vi = get_theta(vi)
    theta2_vi = get_theta2(vi)

    tau_vi = BiCopPar2Tau(family = vi$cop_type, par = theta_vi, par2 = theta2_vi)



    if (vi$structfactor == 1) {
        par(mfrow =c(1,3))
    } else {
        par(mfrow =c(2,3))
    }

    par(mar=c(5,5,3,1))

    hist(v0_vi, xlab = expression(v0[t]))

    hist(tau_vi, xlab = expression(tau[t]))

    hist(theta2_vi, xlab = expression(theta_2[t]))

    if (vi$structfactor == 2 | vi$structfactor == 3){
        latent_theta_vi = get_latent_theta(vi)
        latent_theta2_vi = get_latent_theta2(vi)

        hist(v_vi, xlab = expression(v[t]))

        latent_tau_vi = BiCopPar2Tau(family = vi$latent_copula_type,
                                     par = latent_theta_vi, par2 = latent_theta2_vi)

        hist(latent_tau_vi, xlab = expression(tau_latent[t]) )

        hist(latent_theta2_vi, xlab = expression(theta_latent2[t]))
    }
    if (vi$structfactor == 4) {
        vine_theta_vi = get_vine_theta(vi)
        vine_theta2_vi = get_vine_theta2(vi)
        vine_tau_vi = BiCopPar2Tau(family = vi$vine_copula_type,
                                   par = vine_theta_vi, par2 = vine_theta2_vi)

        hist(vine_tau_vi, xlab = expression(tau_vine[t]) )

        hist(vine_theta2_vi, xlab = expression(theta_vine2[t]))
    }
    par(mfrow =c(1,1))
}

#' @export
plot.hmcfcop <- function(hmc) {
    op <- par('bg', 'col', 'col.axis')
    on.exit(par(op))
    v0_hmc = get_v0(hmc)
    v_hmc = get_v(hmc)


    theta_hmc = get_theta(hmc)
    theta2_hmc = get_theta2(hmc)

    tau_hmc = BiCopPar2Tau(family = hmc$cop_type, par = theta_hmc, par2 = theta2_hmc)

    if (hmc$structfactor == 1) {
        par(mfrow =c(1,3))
    } else {
        par(mfrow =c(2,3))
    }

    par(mar=c(5,5,3,1))

    hist(v0_hmc, xlab = expression(v0[t]))

    hist(tau_hmc, xlab = expression(tau[t]))

    hist(theta2_hmc, xlab = expression(theta_2[t]))

    if (vi$structfactor == 2 | vi$structfactor == 3){
        latent_theta_hmc = get_latent_theta(hmc)
        latent_theta2_hmc = get_latent_theta2(hmc)

        hist(v_hmc, xlab = expression(v[t]))

        vine_tau_hmc = BiCopPar2Tau(family = hmc$vine_copula_type,
                                    par = vine_theta_hmc, par2 = vine_theta2_hmc)

        hist(vine_tau_hmc, xlab = expression(tau_vine[t]) )

        hist(vine_theta2_hmc, xlab = expression(theta_vine2[t]))
    }
    if (hmc$structfactor == 4) {
        vine_theta_hmc = get_vine_theta(hmc)
        vine_theta2_hmc = get_vine_theta2(hmc)

        vine_tau_hmc = BiCopPar2Tau(family = hmc$vine_copula_type,
                                    par = vine_theta_hmc, par2 = vine_theta2_hmc)
        hist(vine_tau_hmc, xlab = expression(tau_vine[t]) )
        hist(vine_theta2_hmc, xlab = expression(theta_vine2[t]))
    }
    par(mfrow =c(1,1))
}


#' @export
compare_vi_hmc <- function(vi,hmc){

    if (! all.equal(vi$structfactor, hmc$structfactor))
        stop("Factor copula model needs to be the same ")

    op <- par('bg', 'col', 'col.axis')
    on.exit(par(op))

    v0_vi = get_v0(vi)
    v_vi = get_v(vi)
    theta_vi = get_theta(vi)
    theta2_vi = get_theta2(vi)

    v0_vi_sd = get_v0_sd(vi)
    theta_vi_sd = get_theta_sd(vi)
    theta2_vi_sd = get_theta2_sd(vi)


    v0_hmc = get_v0(hmc)
    v_hmc = get_v(hmc)
    theta_hmc = get_theta(hmc)
    theta2_hmc = get_theta2(hmc)

    v0_hmc_sd = get_v0_sd(hmc)
    theta_hmc_sd = get_theta_sd(hmc)
    theta2_hmc_sd = get_theta2_sd(hmc)

    # tau_vi = BiCopPar2Tau(family = vi$cop_type, par = theta_vi, par2 = theta2_vi)
    # tau_hmc = BiCopPar2Tau(family = hmc$cop_type, par = theta_hmc, par2 = theta2_hmc)


    if (vi$structfactor == 1) {
        par(mfrow =c(1,4))
    } else {
        par(mfrow =c(2,4))
    }

    par(mar=c(5,5,3,1))

    plot(v0_hmc, v0_vi, xlab = expression(v0[hmc]), ylab = expression(v0[vi]))
    abline(a= 0, b=1, col="red")

    plot(v0_hmc_sd, v0_vi_sd, xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])))
    abline(a= 0, b=1, col="red")

    plot(theta_hmc, theta_vi, xlab = expression(theta[hmc]), ylab = expression(theta[vi]))
    abline(a= 0, b=1, col="red")

    plot(theta_hmc_sd, theta_vi_sd, xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])))
    abline(a= 0, b=1, col="red")

    if (vi$structfactor == 2 | vi$structfactor == 3){
        v_vi_sd = get_v_sd(vi)
        v_hmc_sd = get_v_sd(hmc)

        latent_theta_vi = get_latent_theta(vi)
        latent_theta2_vi = get_latent_theta2(vi)
        latent_theta_vi_sd = get_latent_theta_sd(vi)
        latent_theta2_vi = get_latent_theta2_sd(vi)


        latent_theta_hmc = get_latent_theta(hmc)
        latent_theta2_hmc = get_latent_theta2(hmc)
        latent_theta_hmc_sd = get_latent_theta_sd(hmc)
        latent_theta2_hmc_sd = get_latent_theta2_sd(hmc)

        plot(v_hmc, v_vi, xlab = expression(v[hmc]), ylab = expression(v[vi]))
        abline(a= 0, b=1, col="red")

        plot(v_hmc_sd, v_hmc_sd, xlab = expression(sd(v[hmc]) ), ylab = expression(sd (v[vi])))
        abline(a= 0, b=1, col="red")

        plot(latent_theta_hmc, latent_theta_vi, xlab = expression(theta_latent[hmc]), ylab = expression(theta_latent[vi]))
        abline(a= 0, b=1, col="red")

        plot(latent_theta_hmc_sd, latent_theta_vi_sd, xlab = expression(sd(theta_latent[hmc]) ), ylab = expression(sd (theta_latent[vi])))
        abline(a= 0, b=1, col="red")
    }
    if (vi$structfactor == 4) {

        vine_theta_vi = get_vine_theta(vi)
        vine_theta2_vi = get_vine_theta2(vi)
        vine_theta_vi_sd = get_vine_theta_sd(vi)
        vine_theta2_vi = get_vine_theta2_sd(vi)


        vine_theta_hmc = get_vine_theta(hmc)
        vine_theta2_hmc = get_vine_theta2(hmc)
        vine_theta_hmc_sd = get_vine_theta_sd(hmc)
        vine_theta2_hmc_sd = get_vine_theta2_sd(hmc)

        plot(theta2_hmc, theta2_vi, xlab = expression(theta2[hmc]), ylab = expression(theta2[vi]))
        abline(a= 0, b=1, col="red")

        plot(theta2_hmc_sd, theta2_vi_sd, xlab = expression(sd(theta2[hmc]) ), ylab = expression(sd (theta2[vi])))
        abline(a= 0, b=1, col="red")

        plot(vine_theta_hmc, vine_theta_vi, xlab = expression(theta_vine[hmc]), ylab = expression(theta_vine[vi]))
        abline(a= 0, b=1, col="red")

        plot(vine_theta_hmc_sd, vine_theta_vi_sd, xlab = expression(sd(theta_vine[hmc]) ), ylab = expression(sd (theta_vine[vi])))
        abline(a= 0, b=1, col="red")
    }
    par(mfrow =c(1,1))
}

#' @export
check_edges <- function(vi, vi_rng){
    accurate = 0
    len = nrow(vi$edges) # vi is the accurate structure
    for (i in c(1:len)){
        id = (1:nrow(vi_rng$edges))[(vi_rng$edges[,1] == vi$edges[i,1]) & (vi_rng$edges[,2] == vi$edges[i,2])]
        accurate = accurate + sum( vi_rng$vine_copula_type[id] == vi$vine_copula_type[i] )
    }
    return(accurate)
}
