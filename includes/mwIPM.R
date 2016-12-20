# Setup generic functions
setVars <- function(x) UseMethod("setVars")
setPars <- function(x) UseMethod("setPars")

# Perform fits
setFloweringFit <- function(obj, compute, update) UseMethod("setFloweringFit")
setSurvivalFit <- function(obj, compute, update) UseMethod("setSurvivalFit")
setGrowthFit <- function(obj, compute, update) UseMethod("setGrowthFit")
setPodsFit <- function(obj, compute, update) UseMethod("setPodsFit")

# Set matrices from fits
setFloweringMatrix <- function(x) UseMethod("setFloweringMatrix")
setSurvivalMatrix <- function(x) UseMethod("setSurvivalMatrix")
setHerbivoryMatrix <- function(obj, dist.herb, update) UseMethod("setHerbivoryMatrix")
setPodsMatrix <- function(x) UseMethod("setPodsMatrix")
setGrowthMatrix <- function(x) UseMethod("setGrowthMatrix")
setSeedlingRecruitmentMatrix <- function(x) UseMethod("setSeedlingRecruitmentMatrix")

# Compute kernels from matrices
computeSexualKernel <- function(x) UseMethod("computeSexualKernel")
computeClonalKernel <- function(x) UseMethod("computeClonalKernel")
computeFullKernel <- function(x) UseMethod("computeFullKernel")

# Compute MPM
computeMPM <- function(x) UseMethod("computeMPM")

# Analyze growth rate
analyzeGrowthRate <- function(x) UseMethod("analyzeGrowthRate")

glmerCtrl <- glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun=50000))

# Constructor for mwMod class
mwIPM <- function(x) {
  #
  # Specify min and max of variables in list form.
  #
  if (all(names(x) != "data")) {
    stop("You must specify data (data) for class mwIPM.")
  }
  if (all(names(x) != "N")) {
    x$N = 50
  }
  if (all(names(x) != "site")) {
    x$site = "Bertha"
  }
  if (all(names(x) != "h_apical")) {
    h_apical <- list(min = 0,
                     max = 1.1*max(c(x$data$h_apical),na.rm=T))
  }
  if (all(names(x) != "h_apical.next")) {
    h_apical.next <- list(min = 0,
                          max = 1.3*max(x$data$h_apical.next,na.rm=T))
  }
  if (all(names(x) != "log_herb_avg")) {
    log_herb_avg <- list(min = log(0.1),
                         max = log(6.1))
  }
  
  x <- append(x, list(data_orig = data,
                      vars = list(h_apical = append(h_apical, list(b = NA,
                                                                   x = NA,
                                                                   dx = NA)),
                                  h_apical.next = append(h_apical.next, list(b = NA,
                                                                             x = NA,
                                                                             dx = NA)),
                                  log_herb_avg = append(log_herb_avg, list(b = NA,
                                                                           x = NA,
                                                                           dx = NA))),
                      pars = list(flower.fit = NA,
                                  growth.fit = NA,
                                  surv.fit = NA,
                                  pods.fit = NA,
                                  seedling.fit = NA,
                                  budling.fit = NA,
                                  munched.fit = NA,
                                  budlings.per.stem.fit = NA,
                                  seedling.emergence = NA,
                                  seeds.per.pod = NA,
                                  dist.herb = NA),
                      matrices = list(F = NA,
                                      S = NA,
                                      G = NA,
                                      P = NA,
                                      R = NA),
                      kernels = list(Ks = NA,
                                     Kc = NA,
                                     K = NA),
                      MPM = NA))
  
  y <- structure(x, class = "mwIPM")
  y <- setVars(y)
  y <- setPars(y)
  y <- computeSexualKernel(y)
  y <- computeClonalKernel(y)
  y <- computeFullKernel(y)
  y <- computeMPM(y)
  
  return(y)
}

setVars.mwIPM <- function (obj, update = TRUE) {
  N <- obj$N
  
  attach(obj$vars, warn.conflicts = FALSE)
  
  # h_apical
  h_apical$b = h_apical$min + c(0:N)*(h_apical$max - h_apical$min)/N # boundary points
  h_apical$x = 0.5*(h_apical$b[1:N]+h_apical$b[2:(N+1)])
  h_apical$dx = h_apical$b[2]-h_apical$b[1] # class size
  
  # h_apical.next
  h_apical.next$b = h_apical.next$min + c(0:N)*(h_apical.next$max - h_apical.next$min)/N # boundary points
  h_apical.next$x = 0.5*(h_apical.next$b[1:N]+h_apical.next$b[2:(N+1)])
  h_apical.next$dx = h_apical.next$b[2] - h_apical.next$b[1] # class size
  
  # log_herb_avg
  log_herb_avg$b = log_herb_avg$min + c(0:N)*(log_herb_avg$max - log_herb_avg$min)/N # boundary points
  log_herb_avg$x = 0.5*(log_herb_avg$b[1:N] + log_herb_avg$b[2:(N+1)])
  log_herb_avg$dx = log_herb_avg$b[2] - log_herb_avg$b[1] # class size
  
  obj$vars = list(h_apical = h_apical,
                  h_apical.next = h_apical.next,
                  log_herb_avg = log_herb_avg)
  
  detach(obj$vars)
  
  return(obj)
}

setPars.mwIPM <- function(obj, update = TRUE) {
  #
  # For now, just load them from file here.  In the future, this should be done outside this function.
  # The future is now.  Moving these to proper location.

  attach("../../data/calculated/seedlingFit.RData", warn.conflicts = FALSE)
  attach("../../data/calculated/budlingFit.RData", warn.conflicts = FALSE)
  attach("../../data/calculated/munchedFit.RData", warn.conflicts = FALSE)
  attach("../../data/calculated/budlingsPerStemFit.RData", warn.conflicts = FALSE)
  attach("../../data/calculated/seedlingEmergence.RData", warn.conflicts = FALSE)
  seeds_per_pod_data <- read.csv("../../data/seeddata.csv")

  obj$pars <- list(seedling.fit = seedling.fit,
                   budling.fit = budling.fit,
                   munched.fit = munched.fit,
                   budlings.per.stem.fit = budlings.per.stem.fit,
                   seedling.emergence = seedling.emergence,
                   seeds.per.pod = mean(seeds_per_pod_data$total_seed),
                   dist.herb = NA)
  
  compute = FALSE
  
  obj <- obj %>% setFloweringFit(compute = compute, update = FALSE) %>% 
                 setSurvivalFit(compute = compute, update = FALSE) %>%
                 setGrowthFit(compute = compute, update = FALSE) %>%
                 setPodsFit(compute = compute, update = FALSE)

  detach("file:../../data/calculated/seedlingFit.RData")
  detach("file:../../data/calculated/budlingFit.RData")
  detach("file:../../data/calculated/munchedFit.RData")
  detach("file:../../data/calculated/budlingsPerStemFit.RData")
  detach("file:../../data/calculated/seedlingEmergence.RData")

  if (update) {
    obj <- setFloweringMatrix(obj)
    obj <- setSurvivalMatrix(obj)
    obj <- setGrowthMatrix(obj)
    obj <- setPodsMatrix(obj)
    obj <- setHerbivoryMatrix(obj, update=FALSE)
    obj <- setSeedlingRecruitmentMatrix(obj)
  }
  
  return(obj)
}

# Fits ------------------------------

setFloweringFit.mwIPM <- function(obj, compute = FALSE, update = TRUE) {
  if (!file.exists("../../data/calculated/flowerFit.RData") | (compute)) {
    metadata_usc <- metadata %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        !is.na(fec.flower),
                                        !is.na(surv))
    
    metadata_sc <- metadata_usc %>% mutate_each(funs(as.numeric(scale(.))), 
                                                h_apical, 
                                                log_herb_avg)
    
    cat("Computing flowering fit...")
    flower.mdl <- glmer(fec.flower ~ h_apical + h_apical:log_herb_avg - 1 + (1|site/transect)+(h_apical+log_herb_avg|year), data=metadata_sc, nAGQ=1, family=binomial(), control=glmerCtrl)
    cat("done!\n")
    
    flower.fit <- mwMod(list(mdl = flower.mdl, 
                             vars = c("h_apical", "log_herb_avg"), 
                             scaled = list(h_apical = scale(metadata_usc$h_apical), 
                                           log_herb_avg = scale(metadata_usc$log_herb_avg))))
    
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(flower.fit)
  
    save(flower.fit, file = "../../data/calculated/flowerFit.RData")
  } else {
    load("../../data/calculated/flowerFit.RData")
  }
  
  obj$pars$flower.fit <- flower.fit
  
  return(obj)
}

setSurvivalFit.mwIPM <- function(obj, compute = FALSE, update = TRUE) {
  if (!file.exists("../../data/calculated/survivalFit.RData") | (compute)) {
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        !is.na(surv))
    
    metadata_sc <- metadata_usc %>% mutate_each(funs(as.numeric(scale(.))), 
                                                h_apical, 
                                                log_herb_avg)
    
    cat("Computing survival fit...")
    surv.mdl <- glmer(surv ~ 1 + (log_herb_avg|site/transect)+(1|year), data=metadata_sc, family=binomial(), nAGQ=1, control=glmerCtrl)
    cat("done!\n")
    
    surv.fit <- mwMod(list(mdl = surv.mdl, 
                           vars = c("h_apical", "log_herb_avg"), 
                           scaled = list(h_apical = scale(metadata_usc$h_apical), 
                                         log_herb_avg = scale(metadata_usc$log_herb_avg))))
    
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(surv.fit)
    
    save(surv.fit, file = "../../data/calculated/survivalFit.RData")
  } else {
    load("../../data/calculated/survivalFit.RData")
  }
  
  obj$pars$surv.fit <- surv.fit
  
  return(obj)
}

setGrowthFit.mwIPM <- function(obj, compute = FALSE, update = TRUE) {
  if (!file.exists("../../data/calculated/growthFit.RData") | (compute)) {
    # We only want stems that flowered and survived
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        surv == 1)
    
    metadata_sc <- metadata_usc %>% mutate_each(funs(as.numeric(scale(.))), 
                                                h_apical, 
                                                h_apical.next, 
                                                log_herb_avg)
    
    
    cat("Computing growth fit...")
    growth.mdl <- lmer(h_apical.next ~ h_apical + h_apical:log_herb_avg + (h_apical+log_herb_avg|site/transect)+(h_apical|year), data=metadata_sc, REML=T)
    cat("done!\n")
    
    growth.fit <- mwMod(list(mdl = growth.mdl, 
                             vars = c("h_apical", "log_herb_avg"), 
                             scaled = list(h_apical = scale(metadata_usc$h_apical),
                                           h_apical.next = scale(metadata_usc$h_apical.next),
                                           log_herb_avg = scale(metadata_usc$log_herb_avg))))
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(growth.fit)
    
    save(growth.fit, file = "../../data/calculated/growthFit.RData")
  } else {
    load("../../data/calculated/growthFit.RData")
  }
  
  obj$pars$growth.fit <- growth.fit
  
  return(obj)
}

setPodsFit.mwIPM <- function(obj, compute = FALSE, update = TRUE) {
  if (!file.exists("../../data/calculated/podsFit.RData") | (compute)) {
    # We only want stems that flowered, survived and have data for pods
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        surv == 1,
                                        !is.na(N_pods))
    
    metadata_sc <- metadata_usc %>% mutate_each(funs(as.numeric(scale(.))), 
                                                h_apical, 
                                                h_apical.next, 
                                                log_herb_avg)
    
    cat("Computing pods fit...")
    pods.mdl <- glmer(N_pods ~ h_apical.next + log_herb_avg - 1 + (h_apical.next+log_herb_avg|site/transect)+(h_apical.next+log_herb_avg|year), data=metadata_sc, nAGQ=1, family=poisson(), control=glmerCtrl)
    cat("done!\n")
    
    pods.fit <- mwMod(list(mdl = pods.mdl, 
                           vars = c("h_apical.next", "log_herb_avg"), 
                           scaled = list(h_apical = scale(metadata_usc$h_apical),
                                         h_apical.next = scale(metadata_usc$h_apical.next),
                                         log_herb_avg = scale(metadata_usc$log_herb_avg))))
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(pods.fit)
    
    save(pods.fit, file = "../../data/calculated/podsFit.RData")
  } else {
    load("../../data/calculated/podsFit.RData")
  }
  
  obj$pars$pods.fit <- pods.fit
  
  return(obj)
}

# Matrices ------------------------------

setFloweringMatrix.mwIPM <- function(obj, update = TRUE) {
  # Flowering
  # N x N^2
  obj$matrices$F = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site)})))[rep(1,obj$N),]

  return(obj)
}

setSurvivalMatrix.mwIPM <- function(obj, update = TRUE) {
  # Survival
  # N x N^2
  obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site)})))[rep(1,obj$N),]
  
  return(obj)
}

setGrowthMatrix.mwIPM <- function(obj, update = TRUE) {
  # Growth
  # N x N^2

  N <- obj$N
  Mu <- c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$growth.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site)}))
  G <- matrix(rep(0, N^3), nrow=N)
  for (i in (1:N^2)) {
    for (j in 1:N) {
      G[j,i] = pnorm(obj$vars$h_apical.next$b[j+1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]]) - pnorm(obj$vars$h_apical.next$b[j], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]])
    }
    G[1,i] = G[1,i] + pnorm(obj$vars$h_apical.next$b[1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]])
    G[N,i] = G[N,i] + pnorm(obj$vars$h_apical.next$b[N+1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]], lower.tail=FALSE)
    G[,i] <- G[,i]/obj$vars$h_apical.next$dx
  }
  obj$matrices$G <- G
  
  return(obj)
}

# Check for eviction
# colSums(G%*%H*dx.h_apical.next*dx.log_herb_avg)

# Plot growth
# image.plot(h_apical, h_apical.next, t(G%*%H), col=topo.colors(100))

setPodsMatrix.mwIPM <- function(obj, update = TRUE) {
  # Pods
  # N x N^2

  N <- obj$N
  P <- matrix(rep(0, N^3), nrow=N)
  for (i in 1:N) {
    # h_apical.next
    for (j in 1:N) {
      # log_herb_avg
      P[i,(1+(j-1)*N):(j*N)] <- predict(obj$pars$pods.fit, newdata = data.frame(h_apical.next = obj$vars$h_apical.next$x[i], log_herb_avg = obj$vars$log_herb_avg$x[j]), type=obj$site)
    }
  }
  obj$matrices$P <- P
  
  return(obj)
}

setHerbivoryMatrix.mwIPM <- function(obj, dist.herb = NA, update = TRUE) {
  # Herbivory matrix
  # N x N^2
  
  N <- obj$N
  if (any(is.na(dist.herb))) {
    dist.herb <- obj$pars$munched.fit[[obj$site]]$predict(obj$vars$log_herb_avg$b)
  }
  # Point mass to test growth
  # dist.herb = c(1/dx.log_herb_avg, rep(0,N-1))
  # dist.herb = c(rep(0,N-1), 1/dx.log_herb_avg)
  
  H = matrix(rep(0,N^3), nrow = N^2)
  for (i in 1:N) {
    for (I in 1:N) {
      H[I+N*(i-1), I] = dist.herb[i]
    }
  }

  obj$pars$dist.herb <- dist.herb
  obj$matrices$H <- H
  
  if (update) {
    obj <- computeSexualKernel(obj)
    obj <- computeClonalKernel(obj)
    obj <- computeFullKernel(obj)
  }
  
  return(obj)
}

## Seeding Recruitment
setSeedlingRecruitmentMatrix.mwIPM <- function(obj, update = TRUE) {
  # N x N
  
  R <- t(t(obj$pars$seedling.fit$predict(obj$vars$h_apical$b)))%*%t(rep(1,obj$N))
  
  obj$matrices$R <- R
  
  return(obj)
}

# Kernels ------------------------------

computeSexualKernel.mwIPM <- function(obj) {
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)
  
  Ks <- seedling.emergence[[obj$site]]*seeds.per.pod*R%*%(P*G*S*F)%*%H*log_herb_avg$dx*h_apical.next$dx*h_apical$dx
  # image.plot(h_apical, h_apical, t(Ks), col=topo.colors(100))
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$kernels$Ks <- Ks
  
  return(obj)
}

computeClonalKernel.mwIPM <- function(obj) {
  ## Budling Recruitment (same as Kc)
  
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)
  
  mean.log_herb_avg <- sum(dist.herb*log_herb_avg$x)*log_herb_avg$dx
  mean.buds.per.stem <- exp(predict(budlings.per.stem.fit, new = data.frame(log_herb_mean = mean.log_herb_avg)))
  Kc <- t(t(mean.buds.per.stem*budling.fit[[obj$site]]$predict(h_apical$b)))%*%t(rep(1,obj$N))*h_apical$dx
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$kernels$Kc <- Kc
  
  return(obj)
}

## Full kernel
computeFullKernel.mwIPM <- function(obj) {
  obj$kernels$K <- obj$kernels$Ks + obj$kernels$Kc
  # image.plot(h_apical, h_apical, t(K), col=topo.colors(100))
  # contour(h_apical, h_apical, t(K), add = TRUE, drawlabels = TRUE)
  
  return(obj)
}

# Analysis ------------------------------

computeMPM.mwIPM <- function(obj) {
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)
  
  mean.log_herb_avg <- sum(dist.herb*log_herb_avg$x)*log_herb_avg$dx
  mean.buds.per.stem <- exp(predict(budlings.per.stem.fit, new = data.frame(log_herb_mean = mean.log_herb_avg)))
  
  Kss <- seedling.emergence[[obj$site]]*seeds.per.pod*(P*G*S*F)%*%H*log_herb_avg$dx*h_apical.next$dx*h_apical$dx
  alpha <- sum(Kss%*%t(t(seedling.fit$predict(h_apical$b))))
  beta <- sum(Kss%*%t(t(budling.fit[[obj$site]]$predict(h_apical$b))))
  pems <- mean.buds.per.stem
  pemb <- mean.buds.per.stem
  # pems <- sum(Kc%*%t(t(seedling.fit$predict(h_apical.b))))*dx.h_apical
  # pemb <- sum(Kc%*%t(t(budling.fit[['Bertha']]$predict(h_apical.b))))*dx.h_apical
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$MPM <- matrix(c(alpha, beta, pems, pemb), nrow = 2, byrow=TRUE)
  
  return(obj)
}

analyzeGrowthRate.mwIPM <- function(obj) {
  return(Re(eigen(obj$kernels$K)$values[1]))
}