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
rtheta <-  function(family,tau_gen = 0.9, theta = TRUE) {

    if (theta){
        max_theta = BiCopTau2Par(family, tau_gen, check.taus = TRUE)
        if (family == 1) theta_gen = runif(1, min = 0.2, max = max_theta)
        if (family == 2) theta_gen = runif(1, min = 0.2, max = max_theta)
        if (family == 3) theta_gen = runif(1, min = 0.2, max = max_theta)
        if (family == 4) theta_gen = runif(1, min = 1.1, max = max_theta)
        if (family == 5) theta_gen = runif(1, min = 1, max = max_theta)
        if (family == 6) theta_gen = runif(1, min = 1.2, max = max_theta)

    } else {
        if (family == 1) theta_gen = 0
        if (family == 2) theta_gen = runif(1, min = 2, max = 10)
        if (family == 3) theta_gen = 0
        if (family == 4) theta_gen = 0
        if (family == 5) theta_gen = 0
        if (family == 6) theta_gen = 0
    }
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
fcopsim <- function(t_max, n_max, k_max = 1, family, gid = rep(1,n_max), structfactor = 1, seed_num = 0) {
    set.seed(seed_num)

    if (! all.equal(t_max, as.integer(t_max)))
        stop("'t_max' has to be a integer")
    if (! all.equal(n_max, as.integer(n_max)))
        stop("'n_max' has to be a integer")
    if (! all.equal(k_max, as.integer(k_max)))
        stop("'n_max' has to be a integer")
    if (! (t_max > 0) & (n_max > 0) & (k_max > 0) )
        stop("'t_max', 'n_max', 'k_max' has to be greater than 0")

    if (!(length(family) %in% c(1, n_max)))
        stop("'family' has to be a single number or a size n_max vector")

    if (length(family) == 1){
        family = matrix(family, nrow = n_max, ncol = k_max)
        }

    u <- matrix(runif(t_max*n_max),nrow=t_max, ncol = n_max)
    theta <- matrix(0,nrow=n_max, ncol = k_max)
    theta2 <- matrix(0,nrow=n_max, ncol = k_max)
    tau_mat <- matrix(0,nrow=n_max, ncol = k_max)

    if (structfactor == 1)
    {
        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        for (k in 1:k_max){
            for (i in 1:n_max){
                if (k > 1) {
                    tau_gen = min(tau_mat[i,1:(k-1)])
                } else {
                        tau_gen = 0.65 }

                theta[i,k] <- rtheta(family[i,k],tau_gen, TRUE)
                theta2[i,k] <- rtheta(family[i,k],0, FALSE)
            }
            tau_mat[,k] <- BiCopPar2Tau(family = family[,k], par = theta[,k], par2 = theta2[,k] )
        }

        for (k in k_max:1){
            for (i in 1:n_max){
                obj <- BiCop(family = family[i,k], par = theta[i,k], par2 = theta2[i,k])
                u[,i] <- BiCopHinv2(u[,i], v[,k], obj)
            }
        }


    }

    if (structfactor == 2)
    {
        if (! all.equal(k_max, max(gid)+1))
            stop("'k_max' has to be max(gid)+1")

        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        for (i in 1:n_max){
            # common factor
            tau_gen = 0.65
            theta[i,1] <- rtheta(family[i,1], tau_gen, TRUE)
            theta2[i,1] <- rtheta(family[i,1],0, FALSE)
            tau_mat[i,1] <- BiCopPar2Tau(family = family[i,1], par = theta[i,1], par2 = theta2[i,1] )

            # group factor
            tau_gen = tau_mat[i,1]
            theta[i,2] <- rtheta(family[i,gid[i]+1],tau_gen, TRUE)
            theta2[i,2] <- rtheta(family[i,gid[i]+1],0, FALSE)
            tau_mat[i,2] <- BiCopPar2Tau(family = family[i,gid[i]+1], par = theta[i,2], par2 = theta2[i,2] )

            # group factor
            obj <- BiCop(family = family[i,gid[i]+1], par = theta[i,2], par2 = theta2[i,2])
            u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)

            # common factor
            obj <- BiCop(family = family[i,1], par = theta[i,1], par2 = theta2[i,1])
            u[,i] <- BiCopHinv2(u[,i], v[,1], obj)

        }

    }

    if (structfactor == 3)
    {
        if (! all.equal(k_max, max(gid)+1))
            stop("'k_max' has to be max(gid)+1")

        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        for (k in 2:k_max){
            obj <- BiCop(family = 1, par = 0.7)
            # common factor is v[,1]
            # group factor are v[,i]
            v[,k] <- BiCopHinv2(v[,k], v[,1], obj)
        }

        for (i in 1:n_max){

            tau_gen = 0.65

            theta[i,2] <- rtheta(family[i,gid[i]+1],tau_gen, TRUE)
            theta2[i,2] <- rtheta(family[i,gid[i]+1],0, FALSE)
            tau_mat[i,2] <- BiCopPar2Tau(family = family[i,gid[i]+1], par = theta[i,2], par2 = theta2[i,2] )

            # group factor
            obj <- BiCop(family = family[i,gid[i]+1], par = theta[i,2], par2 = theta2[i,2])
            u[,i] <- BiCopHinv2(u[,i], v[,gid[i]+1], obj)

        }

    }


    list( u = u,
          v = v,
          t_max = t_max,
          n_max = n_max,
          k_max = k_max,
          family = family,
          theta = theta,
          theta2 = theta2,
          gid = gid,
          structfactor = structfactor)
}
