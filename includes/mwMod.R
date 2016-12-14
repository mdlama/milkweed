# Constructor for mwMod class
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

# Method for mwMod prediction
predict.mwMod <- function(obj, newdata, type="Bertha") {
  pars <- obj$pars$unscaled
  linpart <- pars[type,"(Intercept)"] + pars[type,obj$vars[1]]*newdata[[obj$vars[1]]] + pars[type,obj$vars[2]]*newdata[[obj$vars[2]]] + pars[type,paste0(obj$vars[1],":",obj$vars[2])]*newdata[[obj$vars[1]]]*newdata[[obj$vars[2]]]
  if (class(obj$mdl) == "glmerMod") {
    predict.vrate <- obj$mdl@resp$family$linkinv(linpart)
  } else {
    predict.vrate <- linpart
  }
}

calcPars <- function(mdl, scaled, vars) {
  growth <- (class(mdl) == "lmerMod")
  pars <- list(scaled=NA, unscaled=NA)
  pars$scaled <- matrix(rep(0,6*4), byrow=TRUE, nrow=6)
  colnames(pars$scaled) <- c("(Intercept)", vars[1], vars[2], paste0(vars[1], ":", vars[2]))
  rownames(pars$scaled) <- c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")
  
  # Assign values for Bertha
  assignme <- rownames(coef(summary(mdl)))
  for (i in 1:length(assignme)) {
    pars$scaled["Bertha", assignme[i]] <- coef(summary(mdl))[assignme[i], "Estimate"]
  }
  
  # Assign values for sites
  for (i in 1:ncol(coef(mdl)$site)) {
    varexp <- colnames(coef(mdl)$site)[i]
    pars$scaled[2:6,varexp] <- coef(mdl)$site[,varexp]
  }
  
  mur <- 0
  sdr <- 1
  if (growth) {
    mur <- attributes(scaled$h_apical.next)$'scaled:center'
    sdr <- attributes(scaled$h_apical.next)$'scaled:scale'
    
    pars$sd <- rep(0,6)
    names(pars$sd) <- c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")
    pars$sd['Bertha'] <- sdr*sd(mdl@frame$h_apical.next - predict(mdl, type="response", re.form=NA))
    pars$sd[2:6] <- sdr*sd(mdl@frame$h_apical.next - predict(mdl, type="response", re.form=~(h_apical+log_herb_avg|site)))
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
    cbind(rep(mur, 6), matrix(rep(0, 18), nrow=6))
  colnames(pars$unscaled) <- c("(Intercept)", vars[1], vars[2], paste0(vars[1], ":", vars[2]))
  
  calcPars <- pars
}

checkPars <- function(obj) {
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
  sites <- rownames(obj$pars$unscaled)[-1]
  # Get random effect terms
  terms <- sapply(findbars(formula(obj$mdl)), deparse)
  # Get location of site random effect (space is necessary in front)
  I <- which(grepl(" site", terms))
  # Reformulate with ~
  re.form <- reformulate(paste0("(", terms[I], ")"))
  fitresp <- mur + sdr*predict(obj$mdl, type="response", re.form=re.form)
  E <- rep(0, 5)
  for (i in 1:length(sites)) {
    isite <- sites[i]
    I <- which(obj$mdl@frame$site == isite)
    newsitedata <- newdata[I,]
    estresp <- predict(obj, newdata=newsitedata, type=isite)
    E[i] <- max(abs(estresp-fitresp[I]))
  }
  cat(sprintf("Max site errors: %g\n", max(E)))
}