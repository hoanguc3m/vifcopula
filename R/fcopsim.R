#' Sum of vector elements.
#'
#' \code{sum} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }
#' @export
rtheta <-  function(family, tau_min = 0.2, tau_max = 0.8) {
    theta_gen = c(0,0)

    if ((family == 23) | (family == 33) |
            (family == 24) | (family == 34) |
            (family == 26) | (family == 36) ) {
        tau_neg_max = - abs(tau_min)
        tau_neg_min = - abs(tau_max)
        tau_max = tau_neg_max
        tau_min = tau_neg_min
    }

    max_theta = BiCopTau2Par(family, tau_max, check.taus = TRUE)
    min_theta = BiCopTau2Par(family, tau_min, check.taus = TRUE)

    if (family == 1) theta_gen[1] = runif(1, min = min_theta, max = max_theta)
    if (family == 2) {
        unif_rand = BiCopSim(1,family = 1, par = 0.8)
        theta_gen[1] = qunif(unif_rand[1], min = min_theta, max = max_theta)
        theta_gen[2] = qunif(unif_rand[2], min = 2, max = 15)
    }
    if ((family == 3) | (family == 13) | (family == 23) | (family == 33))
        theta_gen[1] = runif(1, min = min_theta, max = max_theta)
    if ((family == 4) | (family == 14) | (family == 24) | (family == 34))
        theta_gen[1] = runif(1, min = min_theta, max = max_theta)
    if (family == 5)
        theta_gen[1] = runif(1, min = min_theta, max = max_theta)
    if ((family == 6) | (family == 16) | (family == 26) | (family == 36))
        theta_gen[1] = runif(1, min = min_theta, max = max_theta)

    theta_gen
}

#' Sum of vector elements.
#'
#' \code{fcopsim} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }
#' @export
fcopsim <- function(t_max, n_max, k_max = 1, family, family_latent = family, gid = rep(1,n_max),
                    tau_range = c(0.2,0.8), tau_latent_range = c(0.2,0.8),
                    structfactor = 1, seed_num = 0) {
    set.seed(seed_num)

    if (! all.equal(t_max, as.integer(t_max)))
        stop("'t_max' has to be a integer")
    if (! all.equal(n_max, as.integer(n_max)))
        stop("'n_max' has to be a integer")
    if (! all.equal(k_max, as.integer(k_max)))
        stop("'n_max' has to be a integer")
    if (! (t_max > 0) & (n_max > 0) & (k_max > 0) )
        stop("'t_max', 'n_max', 'k_max' has to be greater than 0")
    if ( (tau_range[1] > tau_range[2]) )
        stop("Tau correlation range is invalided")
    if ( (tau_latent_range[1] > tau_latent_range[2]) )
        stop("Tau correlation range of the latent is invalided")

    if ( (structfactor == 1) && (k_max > 1))
        stop("Only support one factor model, for two factor model try: structfactor = 2 ")

    # Check copula family
    if (!(length(family) %in% c(1, n_max)))
        stop("'family' has to be a single number or a size n_max vector")
    if (length(family) == 1){
        family = rep(family, n_max)
    } else {
        family = as.vector(family)
    }

    # Check copula latent family
    family_latent = as.vector(family_latent)
    if (structfactor == 1) family_latent = NULL
    if (structfactor == 2) {
        if (length(family_latent) == 1){
            family_latent = rep(family_latent, n_max)
        } else {
            if (length(family_latent) != n_max )
                stop("'family_latent' has to be a single number or a size n_max vector")
        }
        if ( k_max != max(gid)+1)
            stop("'k_max' has to be max(gid)+1")
    }
    if (structfactor == 3) {
        if (length(family_latent) == 1){
            family_latent = rep(family_latent, k_max-1)
        } else {
            if ( length(family_latent) != k_max-1 )
                stop("'family_latent' has to be a single number or a size k_max-1 vector")
        }
        if ( k_max != max(gid)+1)
            stop("'k_max' has to be max(gid)+1")
    }


    u <- matrix(runif(t_max*n_max),nrow=t_max, ncol = n_max)
    theta <- matrix(0,nrow = n_max, ncol = 2)
    theta_latent <- NULL

    tau_mat <- matrix(0,nrow=n_max, ncol = 2)

    if (structfactor == 1)
    {
        v <- matrix(runif(t_max*1), nrow = t_max, ncol = 1)
        for (i in 1:n_max){
            theta[i,] <- rtheta(family[i],tau_range[1], tau_range[2])
            obj <- BiCop(family = family[i], par = theta[i,1], par2 = theta[i,2])
            u[,i] <- BiCopHinv2(u[,i], v[,1], obj)
        }
        tau_mat[,1] <- BiCopPar2Tau(family = family, par = theta[,1], par2 = theta[,2] )
    }

    if (structfactor == 2)
    {

        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        theta_latent <- matrix(0,nrow = n_max, ncol = 2)

        for (i in 1:n_max){
            theta_latent[i,] <- rtheta(family_latent[i],tau_latent_range[1], tau_latent_range[2])

            obj <- BiCop(family = family_latent[i], par = theta_latent[i,1], par2 = theta_latent[i,2])
            # common factor is v[,1]
            # group factor are v[,k+1]
            u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)

        }
        tau_mat[,2] <- BiCopPar2Tau(family = family_latent, par = theta_latent[,1], par2 = theta_latent[,2] )

        for (i in 1:n_max){
            theta[i,] <- rtheta(family[i], tau_range[1], tau_range[2])

            # group factor
            obj <- BiCop(family = family[i], par = theta[i,1], par2 = theta[i,2])
            u[,i] <- BiCopHinv2(u[,i], v[,1], obj)

        }
        tau_mat[,1] <- BiCopPar2Tau(family = family, par = theta[,1], par2 = theta[,2] )

    }

    if (structfactor == 3)
    {
        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        theta_latent <- matrix(0,nrow = k_max-1, ncol = 2)
        for (k in 2:k_max){
            theta_latent[k-1,] <-  rtheta(family_latent[k-1],tau_latent_range[1], tau_latent_range[2])
            obj <- BiCop(family = family_latent[k-1], par = theta_latent[k-1,1], par2 = theta_latent[k-1,2])
            # common factor is v[,1]
            # group factor are v[,k]
            v[,k] <- BiCopHinv2(v[,k], v[,1], obj)
        }
        tau_mat[1:(k_max-1),2] <- BiCopPar2Tau(family = family_latent, par = theta_latent[,1], par2 = theta_latent[,2] )

        for (i in 1:n_max){

            theta[i,] <- rtheta(family[i], tau_range[1], tau_range[2])

            # group factor
            obj <- BiCop(family = family[i], par = theta[i,1], par2 = theta[i,2])
            u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)

        }
        tau_mat[,1] <- BiCopPar2Tau(family = family, par = theta[,1], par2 = theta[,2] )
    }


    list( u = u,
          v = v,
          t_max = t_max,
          n_max = n_max,
          k_max = k_max,
          family = family,
          theta = theta[,1],
          theta2 = theta[,2],
          family_latent = family_latent,
          theta_latent = theta_latent[,1],
          theta2_latent = theta_latent[,2],
          gid = gid,
          structfactor = structfactor)
}


#' #' @export
#' fcoppred <- function(vi, iteration, seed_num = 0) {
#'     set.seed(seed_num)
#'     n_max = vi$n_max
#'     k_max = vi$k_max
#'     family = vi$cop_type
#'     family_latent = vi$latent_copula_type
#'     gid = vi$gid
#'     structfactor = vi$structfactor
#'     theta = get_theta(vi)
#'     theta2 = get_theta2(vi)
#'     theta_latent = get_latent_theta(vi)
#'     theta_latent2 = get_latent_theta2(vi)
#'
#'     if (! all.equal(iteration, as.integer(iteration)))
#'         stop("'iteration' has to be a integer")
#'     if (! all.equal(n_max, as.integer(n_max)))
#'         stop("'n_max' has to be a integer")
#'     if (! all.equal(k_max, as.integer(k_max)))
#'         stop("'n_max' has to be a integer")
#'     if (! (iteration > 0) & (n_max > 0) & (k_max > 0) )
#'         stop("'iteration', 'n_max', 'k_max' has to be greater than 0")
#'
#'     if ( (structfactor == 1) && (k_max > 1))
#'         stop("Only support one factor model, for two factor model try: structfactor = 2 ")
#'
#'
#'     # Check copula latent family
#'     u <- matrix(runif(iteration*n_max),nrow=iteration, ncol = n_max)
#'
#'     if (structfactor == 1)
#'     {
#'         v <- matrix(runif(iteration*1), nrow = iteration, ncol = 1)
#'         for (i in 1:n_max){
#'             obj <- BiCop(family = family[i], par = theta[i], par2 = theta2[i])
#'             u[,i] <- BiCopHinv2(u[,i], v[,1], obj)
#'         }
#'     }
#'
#'     if (structfactor == 2)
#'     {
#'         v <- matrix(runif(iteration*k_max), nrow = iteration, ncol = k_max)
#'
#'         for (i in 1:n_max){
#'             obj <- BiCop(family = family_latent[i], par = theta_latent[i], par2 = theta_latent2[i])
#'             # common factor is v[,1]
#'             # group factor are v[,k+1]
#'             u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)
#'
#'         }
#'
#'         for (i in 1:n_max){
#'             # group factor
#'             obj <- BiCop(family = family[i], par = theta[i], par2 = theta2[i])
#'             u[,i] <- BiCopHinv2(u[,i], v[,1], obj)
#'
#'         }
#'     }
#'
#'     if (structfactor == 3)
#'     {
#'         v <- matrix(runif(iteration*k_max), nrow = iteration, ncol = k_max)
#'
#'         for (k in 2:k_max){
#'             obj <- BiCop(family = family_latent[k-1], par = theta_latent[k-1], par2 = theta_latent2[k-1])
#'             # common factor is v[,1]
#'             # group factor are v[,k]
#'             v[,k] <- BiCopHinv2(v[,k], v[,1], obj)
#'         }
#'
#'         for (i in 1:n_max){
#'             # group factor
#'             obj <- BiCop(family = family[i], par = theta[i], par2 = theta2[i])
#'             u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)
#'         }
#'     }
#'
#'
#'     list( u = u,
#'         v = v,
#'         iteration = iteration,
#'         n_max = n_max,
#'         k_max = k_max,
#'         family = family,
#'         theta = theta,
#'         theta2 = theta2,
#'         family_latent = family_latent,
#'         theta_latent = theta_latent,
#'         theta2_latent = theta_latent2,
#'         gid = gid,
#'         structfactor = structfactor)
#' }


#' @export
fcoppar2tau <- function(vi){
    theta_vi = get_theta(vi)
    theta2_vi = get_theta2(vi)
    tau_vi = BiCopPar2Tau(family = vi$cop_type, par = theta_vi, par2 = theta2_vi)

    latent_tau_vi = NULL


    if (vi$structfactor > 1) {
        latent_theta_vi = get_latent_theta(vi)
        latent_theta2_vi = get_latent_theta2(vi)
        latent_tau_vi = BiCopPar2Tau(family = vi$latent_copula_type,
                                     par = latent_theta_vi, par2 = latent_theta2_vi)

    }

    return(list(tau_vi = tau_vi, latent_tau_vi = latent_tau_vi))
}
