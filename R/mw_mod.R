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
predict.mwMod <- function(obj, newdata, type="Bertha", perturb=rep(0,4)) {
  pars <- obj$pars$unscaled
  linpart <- (pars[type,"(Intercept)"]+perturb[1]) + (pars[type,obj$vars[1]]+perturb[2])*newdata[[obj$vars[1]]] + (pars[type,obj$vars[2]]+perturb[3])*newdata[[obj$vars[2]]] + (pars[type,paste0(obj$vars[1],":",obj$vars[2])]+perturb[4])*newdata[[obj$vars[1]]]*newdata[[obj$vars[2]]]
  if (class(obj$mdl) == "glmerMod") {
    predict.vrate <- obj$mdl@resp$family$linkinv(linpart)
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
  sites <- rownames(coef(mdl)$site)
  years <- rownames(coef(mdl)$year)
  num_sites <- length(sites)
  num_years <- length(years)
  num_pars <- 1+num_sites+num_years # Don't forget Bertha!
  growth <- (class(mdl) == "lmerMod")
  pars <- list(scaled=NA, unscaled=NA)
  pars$scaled <- matrix(rep(0,num_pars*4), byrow=TRUE, nrow=num_pars)
  colnames(pars$scaled) <- c("(Intercept)", vars[1], vars[2], paste0(vars[1], ":", vars[2]))
  rownames(pars$scaled) <- c("Bertha", sites, years)

  # Assign values for Bertha
  assignme <- rownames(coef(summary(mdl)))
  for (i in 1:length(assignme)) {
    pars$scaled["Bertha", assignme[i]] <- coef(summary(mdl))[assignme[i], "Estimate"]
  }

  # Assign values for sites
  for (i in 1:ncol(coef(mdl)$site)) {
    varexp <- colnames(coef(mdl)$site)[i]
    pars$scaled[2:(1+num_sites),varexp] <- coef(mdl)$site[,varexp]
  }

  # Assign values for years
  for (i in 1:ncol(coef(mdl)$year)) {
    varexp <- colnames(coef(mdl)$year)[i]
    pars$scaled[(1+num_sites+1):num_pars,varexp] <- coef(mdl)$year[,varexp]
  }

  mur <- 0
  sdr <- 1
  if (growth) {
    mur <- attributes(scaled$h_apical.next)$'scaled:center'
    sdr <- attributes(scaled$h_apical.next)$'scaled:scale'

    pars$sd <- rep(0,num_pars)
    names(pars$sd) <- rownames(pars$scaled)
    pars$sd['Bertha'] <- sdr*sd(mdl@frame$h_apical.next - predict(mdl, type="response", re.form=NA))
    pars$sd[2:(1+num_sites)] <- sdr*sd(mdl@frame$h_apical.next - predict(mdl, type="response", re.form=~(1|site)))
    pars$sd[(1+num_sites+1):num_pars] <- sdr*sd(mdl@frame$h_apical.next - predict(mdl, type="response", re.form=~(1|year)))
  }
  mux <- c(0,0)
  sdx <- c(0,0)
  for (i in 1:2) {
    mux[i] <- attributes(scaled[[vars[i]]])$'scaled:center'
    sdx[i] <- attributes(scaled[[vars[i]]])$'scaled:scale'
  }

  # Calculate unscaled parameters
  pars$unscaled <- pars$scaled %*% diag(c(sdr, sdr/sdx[1], sdr/sdx[2], sdr/(sdx[1]*sdx[2]))) %*%
    cbind(c(1, -1*mux[1], -1*mux[2], mux[1]*mux[2]),
          c(0,         1,         0,     -1*mux[2]),
          c(0,         0,         1,     -1*mux[1]),
          c(0,         0,         0,            1)) +
    cbind(rep(mur, num_pars), matrix(rep(0, 3*num_pars), nrow=num_pars))
  colnames(pars$unscaled) <- c("(Intercept)", vars[1], vars[2], paste0(vars[1], ":", vars[2]))

  return(pars)
}

#' Check to make sure parameter unscaling worked after call to
#' \code{\link{calcPars}} function.
#'
#' @param obj mwMod model object
#'
#' @export
checkPars <- function(obj) {
  sites <- rownames(coef(obj$mdl)$site)
  years <- rownames(coef(obj$mdl)$year)

  growth <- (class(obj$mdl) == "lmerMod")
  mux <- c(0,0)
  sdx <- c(0,0)
  for (i in 1:2) {
    mux[i] <- attributes(obj$scaled[[obj$vars[i]]])$'scaled:center'
    sdx[i] <- attributes(obj$scaled[[obj$vars[i]]])$'scaled:scale'
  }

  mur <- 0
  sdr <- 1
  if (growth) {
    mur <- attributes(obj$scaled$h_apical.next)$'scaled:center'
    sdr <- attributes(obj$scaled$h_apical.next)$'scaled:scale'
  }
  newdata <- data.frame(row.names = 1:nrow(obj$mdl@frame))
  newdata[[obj$vars[1]]] <- mux[1] + sdx[1]*obj$scaled[[obj$vars[1]]]
  newdata[[obj$vars[2]]] <- mux[2] + sdx[2]*obj$scaled[[obj$vars[2]]]

  # First, check Bertha
  estresp <- predict(obj, newdata=newdata, type='Bertha')
  fitresp <- mur + sdr*predict(obj$mdl, type="response", re.form=NA)
  E <- max(estresp - fitresp)
  cat(sprintf("Bertha error: %g\n", E))

  # Now for sites
  num_sites <- length(sites)
  sites <- rownames(obj$pars$unscaled)[2:(1+num_sites)]
  # Get random effect terms
  terms <- sapply(lme4::findbars(formula(obj$mdl)), deparse)
  # Get location of site random effect (space is necessary in front)
  I <- which(grepl(" site", terms))
  # Reformulate with ~
  re.form <- reformulate(paste0("(", terms[I], ")"))
  fitresp <- mur + sdr*predict(obj$mdl, type="response", re.form=re.form)
  E <- rep(0, length(sites)-1)
  for (i in 1:length(sites)) {
    isite <- sites[i]
    I <- which(obj$mdl@frame$site == isite)
    newsitedata <- newdata[I,]
    estresp <- predict(obj, newdata=newsitedata, type=isite)
    E[i] <- max(abs(estresp-fitresp[I]))
  }
  cat(sprintf("Max site errors: %g\n", max(E)))

  # Now for years
  num_years <- length(years)
  sites <- rownames(obj$pars$unscaled)[(1+num_sites+1):(1+num_sites+num_years)]
  # Get random effect terms
  terms <- sapply(lme4::findbars(formula(obj$mdl)), deparse)
  # Get location of site random effect (space is necessary in front)
  I <- which(grepl(" year", terms))
  # Reformulate with ~
  re.form <- reformulate(paste0("(", terms[I], ")"))
  fitresp <- mur + sdr*predict(obj$mdl, type="response", re.form=re.form)
  E <- rep(0, length(years)-1)
  for (i in 1:length(years)) {
    iyear <- years[i]
    I <- which(obj$mdl@frame$year == iyear)
    newyeardata <- newdata[I,]
    estresp <- predict(obj, newdata=newyeardata, type=iyear)
    E[i] <- max(abs(estresp-fitresp[I]))
  }
  cat(sprintf("Max year errors: %g\n", max(E)))
}
