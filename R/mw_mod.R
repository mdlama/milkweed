# Generics

#' Constructor for mwMod class
#'
#' @param x Argument list
#'
#' @return Object of class mwMod
#' @export
mwMod <- function(x) {
  if (all(names(x) != "mdl")) {
    stop("You must specify a model (mdl) for class mwMod.")
  }
  if (all(names(x) != "vars")) {
    stop("You must specify variables (vars) for class mwMod.")
  }
  if (all(names(x) != "scaled")) {
    stop("You must specify scaled variables (scaled) for class mwMod.")
  }
  x <- append(x, list(pars = calcPars(mdl = x$mdl, scaled = x$scaled, vars = x$vars)))
  structure(x, class = "mwMod")
}

#' Method for mwMod prediction
#'
#' @param obj Object of class mwMod.
#' @param newdata Data frame with new data for prediction.
#' @param type Site
#' @param perturb Perturbation vector for sensitivity analysis.
#'
#' @export
predict.mwMod <- function(obj, newdata, type="Bertha", perturb=rep(0,5)) {
  cat("Hello, World!")
  pars <- obj$pars$unscaled
  linpart <- (pars[type,"(Intercept)"]+perturb[1])
  for (i in 1:length(obj$vars)) {
    linpart <- linpart + (pars[type,obj$vars[i]]+perturb[i+1])*newdata[[obj$vars[i]]]
  }

  if (class(obj$mdl)[1] == "glm") {
    predict.vrate <- obj$mdl$family$linkinv(linpart)
  } else {
    predict.vrate <- linpart
  }
}

#' Store scaled and calculate unscaled parameters for models.
#'
#' @param mdl Model
#' @param scaled List of scaled variables (output of `scale` command)
#' @param vars Vector of names for scaled variables
#'
#' @return List of scaled and unscaled parameters.
calcPars <- function(mdl, scaled, vars) {
  num_sites <- 1
  num_vars <- length(vars)
  pars <- list(scaled=NA, unscaled=NA)
  pars$scaled <- matrix(rep(0,num_sites*(1+num_vars)), byrow=TRUE, nrow=num_sites)
  colnames(pars$scaled) <- c("(Intercept)", vars)
  rownames(pars$scaled) <- c("Bertha")

  pars$unscaled <- pars$scaled
  colnames(pars$unscaled) <- c("(Intercept)", vars)
  rownames(pars$unscaled) <- c("Bertha")

  # Assign values for Bertha
  assignme <- rownames(coef(summary(mdl)))
  for (i in 1:length(assignme)) {
    pars$scaled["Bertha", assignme[i]] <- coef(summary(mdl))[assignme[i], "Estimate"]
  }

  mux <- rep(0,num_vars)
  sdx <- rep(0,num_vars)
  for (i in 1:num_vars) {
    mux[i] <- attributes(scaled[[vars[i]]])$'scaled:center'
    sdx[i] <- attributes(scaled[[vars[i]]])$'scaled:scale'
  }

  # Calculate unscaled parameters
  pars$unscaled[1] <- pars$scaled[1]
  for (i in 1:num_vars) {
    pars$unscaled[1] <- pars$unscaled[1] - pars$scaled[1+i]*mux[i]/sdx[i]
    pars$unscaled[1+i] <- pars$scaled[1+i]/sdx[i]
  }

  return(pars)
}

#' Check to make sure parameter unscaling worked after call to
#' \code{\link{calcPars}} function.
#'
#' @param obj mwMod model object
#'
#' @export
checkPars <- function(obj) {
  num_vars <- length(obj$vars)
  mux <- rep(0,num_vars)
  sdx <- rep(0,num_vars)
  for (i in 1:num_vars) {
    mux[i] <- attributes(obj$scaled[[obj$vars[i]]])$'scaled:center'
    sdx[i] <- attributes(obj$scaled[[obj$vars[i]]])$'scaled:scale'
  }

  newdata <- data.frame(row.names = 1:nrow(obj$mdl$data))
  for (i in 1:num_vars) {
    newdata[[obj$vars[i]]] <- mux[i] + sdx[i]*obj$scaled[[obj$vars[i]]]
  }

  # First, check Bertha
  estresp <- predict(obj, newdata=newdata, type='Bertha')
  fitresp <- predict(obj$mdl, type="response", re.form=NA)
  E <- max(estresp - fitresp)
  cat(sprintf("Bertha error: %g\n", E))
}
