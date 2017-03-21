library(rstan)
library(stringr)
library("inline")

stan_shared <- function(cppfile, outfile, verbose=TRUE) {
    cppcode <- readLines(cppfile)
    ## get model name from the cpp code
    ## Assumes only 1 class in the cpp file
    model.name <- na.omit(str_match(cppcode, "class\\s+(.*?)\\s+:"))[ , 2]
    ## Patch cpp code to be used with Rcpp
    newcppcode <-
        paste("#include <rstan/rstaninc.hpp>",
              paste(cppcode, collapse="\n"),
              rstan:::get_Rcpp_module_def_code(model.name),
              sep = "\n")
    ## Write out new cpp code to a tempfile because
    newcppfile  <- tempfile(fileext=".cpp")
    writeLines(newcppcode, newcppfile)
    ## Before compiling Set all the environment variables for linking,
    ## etc.  This is what inline:::cxxfunction appears to do
    settings <- getPlugin("rstan")

    ### Copied from Rcpp::cxxfunction
    if (!is.null(env <- settings$env)) {
        do.call(Sys.setenv, env)
        if (isTRUE(verbose)) {
            cat(" >> setting environment variables: \n")
            writeLines(sprintf("%s = %s", names(env), env))
        }
    }
    LinkingTo <- settings$LinkingTo
    if (!is.null(LinkingTo)) {
        paths <- .find.package(LinkingTo, quiet = TRUE)
        if (length(paths)) {
            flag <- paste(paste0("-I\"", paths, "/include\""),
                collapse = " ")
            Sys.setenv(CLINK_CPPFLAGS = flag)
            if (isTRUE(verbose)) {
                cat(sprintf("\n >> LinkingTo : %s\n", paste(LinkingTo,
                  collapse = ", ")))
                cat("CLINK_CPPFLAGS = ", flag, "\n\n")
            }
        }
    }
    ### End of cxxfunction copied

    ## Also see inline:::compileCode for some platform indep hacks.
    ## Use R CMD SHLIB to compile
    ## Could also do R CMD COMPILE and then R CMD SHLIB
    cmd <- sprintf("-o %s %s", shQuote(outfile), shQuote(newcppfile))
    R <- file.path(R.home(component = "bin"), "R")
    system2(R, c("CMD", "SHLIB",
                 sprintf("-o %s", shQuote(outfile)),
                 shQuote(newcppfile)))
}

## Compile foo.cpp into a shared object foo.so
## This will need to be generalized for cross platform
stan_shared("outfile.cpp", "foo.so")
## Load the so
dyn.load('foo.so')
## foo.so should appear
getLoadedDLLs()

## See section 3.3 of "Exposing C++ functions and classes with Rcpp modules"
## Load Rcpp Module
mod <- Module("foo", getDynLib('foo'))
foo <- `$`(mod, "foo")
stanmodel_object <- new(foo)

##' List classes within a module
##' see getMethods("show", "Module") from which this was extracted
##' as far as I know, there is no public facing code to view the classes / functions in
##' a Module
module_classes <- function(object) {
    pointer <- Rcpp:::.getModulePointer(object, FALSE)
    if (identical(pointer, Rcpp:::.badModulePointer)) {
        object <- as.environment(object)
        txt <- sprintf("Uninitialized module named \"%s\" from package \"%s\"",
            get("moduleName", envir = object), get("packageName",
                envir = object))
        writeLines(txt)
    } else {
        info <- .Call(Rcpp:::Module__classes_info, pointer)
        names(info)
    }
}
