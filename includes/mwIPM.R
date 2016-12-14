# Setup generic functions
setVars <- function(x) UseMethod("setVars")
setPars <- function(x) UseMethod("setPars")

# Set matrices from fits
setFlowering <- function(x) UseMethod("setFlowering")
setSurvival <- function(x) UseMethod("setSurvival")
setHerbivory <- function(obj, dist.herb, update) UseMethod("setHerbivory")
setPods <- function(x) UseMethod("setPods")
setGrowth <- function(x) UseMethod("setGrowth")
setSeedlingRecruitment <- function(x) UseMethod("setSeedlingRecruitment")

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
  #

  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/vitalFits.RData", warn.conflicts = FALSE)
  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/seedlingFit.RData", warn.conflicts = FALSE)
  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/budlingFit.RData", warn.conflicts = FALSE)
  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/munchedFit.RData", warn.conflicts = FALSE)
  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/budlingsPerStemFit.RData", warn.conflicts = FALSE)
  attach("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/seedlingEmergence.RData", warn.conflicts = FALSE)
  seeds_per_pod_data <- read.csv("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/Milkweed/DATA SETS/seeds_per_pod_final.csv")

  obj$pars <- list(flower.fit = flower.fit,
                   growth.fit = growth.fit,
                   surv.fit = surv.fit,
                   pods.fit = pods.fit,
                   seedling.fit = seedling.fit,
                   budling.fit = budling.fit,
                   munched.fit = munched.fit,
                   budlings.per.stem.fit = budlings.per.stem.fit,
                   seedling.emergence = seedling.emergence,
                   seeds.per.pod = mean(seeds_per_pod_data$total_seed),
                   dist.herb = NA)
 
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/vitalFits.RData")
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/seedlingFit.RData")
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/budlingFit.RData")
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/munchedFit.RData")
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/budlingsPerStemFit.RData")
  detach("file:~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/seedlingEmergence.RData")

  if (update) {
    obj <- setFlowering(obj)
    obj <- setSurvival(obj)
    obj <- setGrowth(obj)
    obj <- setPods(obj)
    obj <- setHerbivory(obj, update=FALSE)
    obj <- setSeedlingRecruitment(obj)
  }
  
  return(obj)
}

setFlowering.mwIPM <- function(obj, update = TRUE) {
  # Flowering
  # N x N^2
  obj$matrices$F = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site)})))[rep(1,obj$N),]

  return(obj)
}

setSurvival.mwIPM <- function(obj, update = TRUE) {
  # Survival
  # N x N^2
  obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site)})))[rep(1,obj$N),]
  
  return(obj)
}

setGrowth.mwIPM <- function(obj, update = TRUE) {
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

setPods.mwIPM <- function(obj, update = TRUE) {
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

setHerbivory.mwIPM <- function(obj, dist.herb = NA, update = TRUE) {
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
setSeedlingRecruitment.mwIPM <- function(obj, update = TRUE) {
  # N x N
  
  R <- t(t(obj$pars$seedling.fit$predict(obj$vars$h_apical$b)))%*%t(rep(1,obj$N))
  
  obj$matrices$R <- R
  
  return(obj)
}

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

## MPM
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