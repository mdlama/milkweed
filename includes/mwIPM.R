# Declarations ----------------------

# Setup generic functions
setVars <- function(obj) UseMethod("setVars")
setPars <- function(obj, compute, saveresults, update) UseMethod("setPars")

# Constants
setSeedsPerPodConst <- function(obj, compute, saveresults, update) UseMethod("setSeedsPerPodConst")
setSeedlingEmergenceConst <- function(obj, compute, saveresults, update) UseMethod("setSeedlingEmergenceConst")

# Perform fits

## Vital rates
setFloweringFit <- function(obj, compute, saveresults, update) UseMethod("setFloweringFit")
setSurvivalFit <- function(obj, compute, saveresults, update) UseMethod("setSurvivalFit")
setGrowthFit <- function(obj, compute, saveresults, update) UseMethod("setGrowthFit")
setPodsFit <- function(obj, compute, saveresults, update) UseMethod("setPodsFit")

## Distributions
setSeedlingDistFit <- function(obj, compute, saveresults, update) UseMethod("setSeedlingDistFit")
setBudlingDistFit <- function(obj, compute, saveresults, update) UseMethod("setBudlingDistFit")
setHerbivoryDistFit <- function(obj, compute, saveresults, update) UseMethod("setHerbivoryDistFit")

## Regressions
setBudlingsPerStemFit <- function(obj, compute, saveresults, update) UseMethod("setBudlingsPerStemFit")

# Set matrices from fits
setFloweringMatrix <- function(obj, update, perturb) UseMethod("setFloweringMatrix")
setSurvivalMatrix <- function(obj, update, perturb) UseMethod("setSurvivalMatrix")
setPodsMatrix <- function(obj, update, perturb) UseMethod("setPodsMatrix")
setGrowthMatrix <- function(obj, update, perturb) UseMethod("setGrowthMatrix")
setHerbivoryMatrix <- function(obj, dist.herb, update) UseMethod("setHerbivoryMatrix")
setSeedlingRecruitmentMatrix <- function(obj, update) UseMethod("setSeedlingRecruitmentMatrix")

# Compute kernels from matrices
computeSexualKernel <- function(obj) UseMethod("computeSexualKernel")
computeClonalKernel <- function(obj) UseMethod("computeClonalKernel")
computeFullKernel <- function(obj) UseMethod("computeFullKernel")

# Setup or compute MPM/IPM
computeMPM <- function(obj) UseMethod("computeMPM")
setSite <- function(obj, site, compute) UseMethod("setSite")
bootIPM <- function(obj) UseMethod("bootIPM")

# Analyze growth rate
analyzeGrowthRate <- function(x) UseMethod("analyzeGrowthRate")

# Renderers
renderFloweringFit <- function(obj) UseMethod("renderFloweringFit") 
renderBudlingDistFit <- function(obj) UseMethod("renderBudlingDistFit")
renderHerbivoryDistFit <- function(obj) UseMethod("renderHerbivoryDistFit")

glmerCtrl <- glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun=50000))

# Constructor ----------------------

mwIPM <- function(x = list()) {
  #
  # Specify min and max of variables in list form.
  #
  
  if (all(names(x) != "data")) {
    x$data <- read.csv(mwROOT("data","stemdata.csv"))
  }
  if (all(names(x) != "N")) {
    x$N = 50
  }
  if (all(names(x) != "site")) {
    site = "Bertha"
    x$site <- NA
  } else if (!(site %in% c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB"))) {
    stop(paste(site, "is not a valid site!  Please choose one of Bertha, BLD1, BLD2, PWR, SKY, or YTB."))
  }
  if (all(names(x) != "compute")) {
    compute = FALSE
  } else if (!is.logical(x$compute)) {
    stop(paste("Variable 'compute' must be TRUE or FALSE."))
  } else {
    compute <- x$compute
  }
  if (all(names(x) != "saveresults")) {
    saveresults = FALSE
  } else if (!is.logical(x$saveresults)) {
    stop(paste("Variable 'saveresults' must be TRUE or FALSE."))
  } else {
    saveresults <- x$saveresults
  }
  x <- x[which(!(names(x) %in% c("compute","saveresults")))]
  
  x <- append(x, list(data_orig = x$data,
                      vars = list(h_apical = list(min = NA,
                                                  max = NA,
                                                  b = NA,
                                                  x = NA,
                                                  dx = NA),
                                  h_apical.next = list(min = NA,
                                                       max = NA,
                                                       b = NA,
                                                       x = NA,
                                                       dx = NA),
                                  log_herb_avg = list(min = NA,
                                                      max = NA,
                                                      b = NA,
                                                      x = NA,
                                                      dx = NA)),
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
  
  y <- structure(x, class = "mwIPM") %>% 
       setPars(compute = compute, saveresults = saveresults, update = FALSE) %>% 
       setSite(compute = compute, site = site) %>%
       computeMPM()
  
  return(y)
}

bootIPM.mwIPM <- function (obj) {
  obj$data <- obj$data_orig %>% group_by(site) %>% 
                                sample_frac(replace=TRUE) %>% 
                                ungroup()
  obj %<>% setPars(compute = TRUE, update = FALSE) %>% 
           setSite(compute = TRUE) %>%
           computeMPM()

  return(obj)
}

# Vars & Pars ----------------------

setVars.mwIPM <- function (obj) {
  N <- obj$N
  
  if (obj$site != "Bertha") {
    data <- obj$data %>% filter(site == obj$site)
  }
  data <- obj$data %>% summarize(max_h_apical = max(h_apical, na.rm=T),
                             max_h_apical.next = max(h_apical.next, na.rm=T))
  
  # h_apical
  h_apical <- list(min = 0,
                   max = 1.1*data$max_h_apical)
  h_apical$b = h_apical$min + c(0:N)*(h_apical$max - h_apical$min)/N # boundary points
  h_apical$x = 0.5*(h_apical$b[1:N]+h_apical$b[2:(N+1)])
  h_apical$dx = h_apical$b[2]-h_apical$b[1] # class size
  
  # h_apical.next
  h_apical.next <- list(min = 0,
                        max = 1.3*data$max_h_apical.next)
  h_apical.next$b = h_apical.next$min + c(0:N)*(h_apical.next$max - h_apical.next$min)/N # boundary points
  h_apical.next$x = 0.5*(h_apical.next$b[1:N]+h_apical.next$b[2:(N+1)])
  h_apical.next$dx = h_apical.next$b[2] - h_apical.next$b[1] # class size
  
  # log_herb_avg
  log_herb_avg <- list(min = log(0.1),
                       max = log(6.1))
  log_herb_avg$b = log_herb_avg$min + c(0:N)*(log_herb_avg$max - log_herb_avg$min)/N # boundary points
  log_herb_avg$x = 0.5*(log_herb_avg$b[1:N] + log_herb_avg$b[2:(N+1)])
  log_herb_avg$dx = log_herb_avg$b[2] - log_herb_avg$b[1] # class size
  
  obj$vars = list(h_apical = h_apical,
                  h_apical.next = h_apical.next,
                  log_herb_avg = log_herb_avg)
  
  return(obj)
}

setPars.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  #
  # For now, just load them from file here.  In the future, this should be done outside this function.
  # The future is now.  Moving these to proper location.

  obj$pars <- list(dist.herb = NA)
  
  obj <- obj %>% setSeedsPerPodConst(compute = compute, 
                                     saveresults = saveresults, 
                                     update = FALSE) %>%
                 setSeedlingEmergenceConst(compute = compute, 
                                           saveresults = saveresults, 
                                           update = FALSE) %>%
                 setBudlingsPerStemFit(compute = compute, 
                                       saveresults = saveresults, 
                                       update = FALSE) %>%
                 setFloweringFit(compute = compute, 
                                 saveresults = saveresults, 
                                 update = FALSE) %>% 
                 setSurvivalFit(compute = compute, 
                                saveresults = saveresults, 
                                update = FALSE) %>%
                 setGrowthFit(compute = compute, 
                              saveresults = saveresults, 
                              update = FALSE) %>%
                 setPodsFit(compute = compute, 
                            saveresults = saveresults, 
                            update = FALSE) %>%
                 setSeedlingDistFit(compute = compute, 
                                    saveresults = saveresults, 
                                    update = FALSE) %>%
                 setBudlingDistFit(compute = compute, 
                                   saveresults = saveresults, 
                                   update = FALSE) %>%
                 setHerbivoryDistFit(compute = compute, 
                                     saveresults = saveresults, 
                                     update = FALSE)

  if (update) {
    obj <- obj %>% setFloweringMatrix() %>%
                   setSurvivalMatrix() %>%
                   setGrowthMatrix() %>%
                   setPodsMatrix() %>%
                   setHerbivoryMatrix(update=FALSE) %>%
                   setSeedlingRecruitmentMatrix()
  }
  
  return(obj)
}

# Constants ------------------------------

setSeedsPerPodConst.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","seedsPerPodConst.RData")) | (compute)) {
    seeds_per_pod_data <- read.csv(mwROOT("data","seeddata.csv"))
    
    cat("Computing seeds per pod constant...")
    seeds.per.pod = mean(seeds_per_pod_data$total_seed)
    cat("done!\n")
    
    if (saveresults) {
      save(seeds.per.pod, file = mwROOT("data","calculated","seedsPerPodConst.RData"))
    }
  } else {
    load(mwROOT("data","calculated","seedsPerPodConst.RData"))
  }
  
  obj$pars$seeds.per.pod <- seeds.per.pod
    
  return(obj)
}

setSeedlingEmergenceConst.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","seedlingEmergenceConst.RData")) | (compute)) {  
    sites <- c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")
    seedling.emergence <- rep(NA, 6)
    names(seedling.emergence) <- sites
    
    # Need seeds per pod, so compute if needed
    if (is.na(obj$pars$seeds.per.pod)) {
      cat("Need seeds per pod to compute seedling emergence.  Loading now...")
      obj <- setSeedsPerPod(obj)
      cat("done!\n")
    }
    seeds.per.pod <- obj$pars$seeds.per.pod
    
    # Bertha
    cat("Computing seedling emergence constants...")
    data_gp <- obj$data %>% group_by(year) %>% summarize(N_seedlings = sum(seedling, na.rm=T), 
                                                         N_seeds = seeds.per.pod*sum(N_pods, na.rm=T))
    
    seedling.emergence[1] <- mean(data_gp$N_seedlings[2:3]/data_gp$N_seeds[1:2])
    
    # Sites
    data_gp <- obj$data %>% group_by(year, site) %>% 
                            summarize(N_seedlings = sum(seedling, na.rm=T), 
                                      N_seeds = seeds.per.pod*sum(N_pods, na.rm=T))
    
    data13_14 <- data_gp %>% filter(year %in% 2013:2014) %>% 
                             group_by(site) %>% 
                             summarize(emergence = last(N_seedlings)/first(N_seeds))

    data14_15 <- data_gp %>% filter(year %in% 2014:2015) %>% 
                             group_by(site) %>% 
                             summarize(emergence = last(N_seedlings)/first(N_seeds))
    
    fulldat <- bind_rows(data13_14, data14_15)
    
    seedling.emergence[2:6] <- (fulldat %>% group_by(site) %>% summarize(emergence = mean(emergence)))$emergence
    cat("done!\n")
    
    if (saveresults) {
      save(seedling.emergence, file = mwROOT("data","calculated","seedlingEmergenceConst.RData"))
    }
  } else {
    load(mwROOT("data","calculated","seedlingEmergenceConst.RData"))
  }
  
  obj$pars$seedling.emergence <- seedling.emergence
  
  return(obj)  
}

# Fits ------------------------------

## Vital rates
setFloweringFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","flowerFit.RData")) | (compute)) {
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
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
  
    if (saveresults) {
      save(flower.fit, file = mwROOT("data","calculated","flowerFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","flowerFit.RData"))
  }
  
  obj$pars$flower.fit <- flower.fit
  
  return(obj)
}

setSurvivalFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","survivalFit.RData")) | (compute)) {
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
    
    if (saveresults) {
      save(surv.fit, file = mwROOT("data","calculated","survivalFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","survivalFit.RData"))
  }
  
  obj$pars$surv.fit <- surv.fit
  
  return(obj)
}

setGrowthFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","growthFit.RData")) | (compute)) {
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
    
    if (saveresults) {
      save(growth.fit, file = mwROOT("data","calculated","growthFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","growthFit.RData"))
  }
  
  obj$pars$growth.fit <- growth.fit
  
  return(obj)
}

setPodsFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","podsFit.RData")) | (compute)) {
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
    
    if (saveresults) {
      save(pods.fit, file = mwROOT("data","calculated","podsFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","podsFit.RData"))
  }
  
  obj$pars$pods.fit <- pods.fit
  
  return(obj)
}

## Distributions
setSeedlingDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","seedlingDistFit.RData")) | (compute)) {
    h_apical <- (obj$data %>% filter(seedling == 1))$h_apical
    
    cat("Computing seedling distribution fit...")
    f1 <- fitdist(h_apical, "lnorm")
    cat("done!\n")
    
    seedling.fit <- vector("list", 2)
    seedling.fit[[1]] <- f1
    seedling.fit[[2]] <- eval(parse(
      text = sprintf(
        "function(x) {
           N <- length(x)
           dx <- x[2]-x[1]
           y <- rep(0, N-1)
           for (j in 1:(N-1)) {
             y[j] = p%s(x[j+1], %g, %g) - p%s(x[j], %g, %g)
           }
           y <- y/(dx*sum(y))
         }",
         seedling.fit[[1]]$distname,
         seedling.fit[[1]]$estimate[1],
         seedling.fit[[1]]$estimate[2],
         seedling.fit[[1]]$distname,
         seedling.fit[[1]]$estimate[1],
         seedling.fit[[1]]$estimate[2]
      )
    ))
    names(seedling.fit) <- c("fit", "predict")
    
    if (saveresults) {
      save(seedling.fit, file = mwROOT("data","calculated","seedlingDistFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","seedlingDistFit.RData"))
  }
  
  obj$pars$seedling.fit <- seedling.fit
  
  return(obj)
}

setBudlingDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","budlingDistFit.RData")) | (compute)) {
    sites <- c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")
    mdls <- c("norm", "weibull", "norm", "weibull", "gamma", "lnorm")
    budling.fit <- vector('list', 6)
    names(budling.fit) <- sites
    
    for (i in 1:6) {
      if (i == 1) { # Bertha
        h_apical <- (obj$data %>% filter(seedling == 0,
                                         !is.na(h_apical)))$h_apical
      } else { # Sites
        h_apical <- (obj$data %>% filter(seedling == 0, 
                                         site == sites[i],
                                         !is.na(h_apical)))$h_apical
      }

      cat("Computing budling distribution fit for", sites[i], "...")
      f0 <- fitdist(h_apical, mdls[i])
      cat("done!\n")

      budling.fit[[i]] <- vector("list", 2)
      budling.fit[[i]][[1]] <- f0
      budling.fit[[i]][[2]] <-
        eval(parse(
          text = sprintf(
            "function(x) {
               N <- length(x)
               dx <- x[2]-x[1]
               y <- rep(0, N-1)
               for (j in 1:(N-1)) {
                 y[j] = p%s(x[j+1], %g, %g) - p%s(x[j], %g, %g)
               }
               y <- y/(dx*sum(y))
            }",
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$estimate[1],
            budling.fit[[i]][[1]]$estimate[2],
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$estimate[1],
            budling.fit[[i]][[1]]$estimate[2]
          )
        ))
      names(budling.fit[[i]]) <- c("fit", "predict")
    }
    
    if (saveresults) {
      save(budling.fit, file = mwROOT("data","calculated","budlingDistFit.RData"))
    }
} else {
  load(mwROOT("data","calculated","budlingDistFit.RData"))
}
  
  obj$pars$budling.fit <- budling.fit
  
  return(obj)
}

setHerbivoryDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","herbivoryDistFit.RData")) | (compute)) {
    sites <- c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")
    mdls <- c("lnorm", "lnorm", "gamma", "lnorm", "lnorm", "gamma")
    munched.fit <- vector('list', 6)
    names(munched.fit) <- sites
    
    for (i in 1:6) {
      if (i == 1) { # Bertha
        thissite <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(munched))
      } else { # Sites
        thissite <- obj$data %>% filter(site == sites[i],
                                        !is.na(h_apical),
                                        !is.na(munched))
      }
      pmunch <- sum(thissite$munched == 1)/nrow(thissite)
      herb_avg <- (thissite %>% filter(munched == 1))$herb_avg
      
      cat("Computing herbivory distribution fit for", sites[i], "...")
      f0 <- fitdist(herb_avg, mdls[i])
      cat("done!\n")
      
      munched.fit[[i]] <- vector("list", 3)
      munched.fit[[i]][[1]] <- f0
      munched.fit[[i]][[2]] <- pmunch
      munched.fit[[i]][[3]] <-
        eval(parse(
          text = sprintf(
            "function(x, justmunch=FALSE) {
               N <- length(x)
               dx <- x[2]-x[1]
               z <- exp(x)-0.1
               y <- rep(0, N)
               for (j in 1:(N-1)) {
                 y[j] = p%s(z[j+1], %g, %g) - p%s(z[j], %g, %g)
               }
               y <- (%g/(dx*sum(y)))*y
               if (!justmunch) {
                 y[1] <- y[1] + (1-%g)/dx
               }
               y <- y[1:(N-1)]
            }",
            munched.fit[[i]][[1]]$distname,
            munched.fit[[i]][[1]]$estimate[1],
            munched.fit[[i]][[1]]$estimate[2],
            munched.fit[[i]][[1]]$distname,
            munched.fit[[i]][[1]]$estimate[1],
            munched.fit[[i]][[1]]$estimate[2],
            munched.fit[[i]][[2]],
            munched.fit[[i]][[2]]
          )
        ))
      
      names(munched.fit[[i]]) <- c("fit", "pmunch", "predict")
    }
    
    if (saveresults) {
      save(munched.fit, file = mwROOT("data","calculated","herbivoryDistFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","herbivoryDistFit.RData"))
  }
  
  obj$pars$munched.fit <- munched.fit
  
  return(obj)
}

setBudlingsPerStemFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(mwROOT("data","calculated","budlingsPerStemFit.RData")) | (compute)) {  
    data_gp <- obj$data %>% group_by(year, transect) %>% 
                            summarize(N_seedlings = sum(seedling, na.rm=T), 
                                      N_total = sum(aliveJune, na.rm=T), 
                                      N_budlings = N_total - N_seedlings,
                                      log_herb_mean = mean(log_herb_avg, na.rm=T),
                                      site = first(site))
    
    data13_14 <- data_gp %>% filter(year %in% 2013:2014) %>% 
                             group_by(transect) %>% 
                             summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                                       log_herb_mean = first(log_herb_mean),
                                       site = first(site))
    
    data14_15 <- data_gp %>% filter(year %in% 2014:2015 & transect != 70) %>%
                             group_by(transect) %>% 
                             summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                                       log_herb_mean = first(log_herb_mean),
                                       site = first(site))
    
    fulldat <- bind_rows(data13_14, data14_15)

    merged <- fulldat %>% group_by(transect) %>% 
                          summarize(log_herb_mean = mean(log_herb_mean),
                                    bdlgs_per_stem = mean(bdlgs_per_stem),
                                    site = first(site))
    
    cat("Calculating budlings per stem fit...")
    exp.model <- lm(log(bdlgs_per_stem) ~ log_herb_mean, data=merged)
    cat("done!\n")

    budlings.per.stem.fit <- exp.model
    budlings.per.stem.fit$site <- merged$site
    
    if (saveresults) {
      save(budlings.per.stem.fit, file = mwROOT("data","calculated","budlingsPerStemFit.RData"))
    }
  } else {
    load(mwROOT("data","calculated","budlingsPerStemFit.RData"))
  }
  
  obj$pars$budlings.per.stem.fit <- budlings.per.stem.fit
  
  return(obj)
}

# Matrices ------------------------------

setFloweringMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Flowering
  # N x N^2
  obj$matrices$F = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]

  return(obj)
}

setSurvivalMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Survival
  # N x N^2
  obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]
  
  return(obj)
}

setGrowthMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Growth
  # N x N^2

  N <- obj$N
  Mu <- c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$growth.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site, perturb=perturb)}))
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
  
  # Check for eviction
  # colSums(G%*%H*dx.h_apical.next*dx.log_herb_avg)
  
  # Plot growth
  # image.plot(h_apical, h_apical.next, t(G%*%H), col=topo.colors(100))
  
  return(obj)
}

setPodsMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Pods
  # N x N^2

  N <- obj$N
  P <- matrix(rep(0, N^3), nrow=N)
  for (i in 1:N) {
    # h_apical.next
    for (j in 1:N) {
      # log_herb_avg
      P[i,(1+(j-1)*N):(j*N)] <- predict(obj$pars$pods.fit, newdata = data.frame(h_apical.next = obj$vars$h_apical.next$x[i], log_herb_avg = obj$vars$log_herb_avg$x[j]), type=obj$site, perturb=perturb)
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

computeFullKernel.mwIPM <- function(obj) {
  obj$kernels$K <- obj$kernels$Ks + obj$kernels$Kc
  # image.plot(h_apical, h_apical, t(K), col=topo.colors(100))
  # contour(h_apical, h_apical, t(K), add = TRUE, drawlabels = TRUE)
  
  return(obj)
}

# Analysis ------------------------------

setSite.mwIPM <- function(obj, site = "Bertha", compute = FALSE) {
  if ((compute) | is.na(obj$site) | (site != obj$site)) {
    if (site %in% c("Bertha", "BLD1", "BLD2", "PWR", "SKY", "YTB")) {
      obj$site <- site
      
      obj <- setVars(obj)
      
      obj <- obj %>% 
        setFloweringMatrix() %>% 
        setSurvivalMatrix() %>% 
        setGrowthMatrix() %>% 
        setPodsMatrix() %>%
        setHerbivoryMatrix(update=FALSE) %>% 
        setSeedlingRecruitmentMatrix()
      
      obj <- obj %>%
        computeSexualKernel() %>% 
        computeClonalKernel() %>%
        computeFullKernel()
      
      obj <- computeMPM(obj)
    } else {
      stop(paste(site, "is not a valid site!  Please choose one of Bertha, BLD1, BLD2, PWR, SKY, or YTB."))
    }
  } else {
    warning(paste(site, "is already the current site."))
  }
  
  return(obj)
}

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

# Renderers ------------------------------

renderFloweringFit.mwIPM <- function(obj) {
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$vars, warn.conflicts = FALSE)
  
  ## Top (heatmap)
  M <- mesh(h_apical$x, log_herb_avg$x)
  plotdata <- tbl_df(data.frame(h_apical = as.vector(M$x),
                                log_herb_avg = as.vector(M$y)))
  z <- predict(flower.fit, 
               newdata=data.frame(h_apical = as.vector(M$x),
                                  log_herb_avg = as.vector(M$y)),
               type="Bertha")
  plotdata <- plotdata %>% mutate(prob.flower = z) 
  
  # Heatmap
  imgt <- ggplot(plotdata, aes(x = h_apical, y = log_herb_avg, z = prob.flower)) +
    geom_raster(aes(fill = prob.flower)) +
    geom_contour(colour = "white", alpha = 0.8) + 
    scale_fill_gradientn("Flowering\nProbability", 
                         colours=c("#00000000","#BBBBBBBB"),
                         limits=c(0, 1))
  
  imgt <- imgt + theme(axis.line = element_blank()) + 
    scale_x_continuous(limits = c(0, h_apical$max), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(log_herb_avg$min, log_herb_avg$max), expand = c(0, 0)) +
    ylab("Log(Herbivory)")
  
  # Add lines
  herb_ex <- data.frame(yintercept = c(log(0.1), mean(obj$data$log_herb_avg, na.rm=T), log(6.1)))
  
  imgt <- imgt + geom_hline(aes(yintercept = yintercept),
                          data = herb_ex,
                          linetype = c(1, 2, 5),
                          col = brewer.pal(5, "YlOrRd")[2:4],
                          size=1)
  
  # Move over a bit to match panel (b) below
  imgt <- imgt + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0,
                                    margin=margin(b = 15, unit = "pt")))

  ## Bottom (curves)
  min <- data.frame(h_apical = h_apical$x, 
                    prob.flower = predict(flower.fit, 
                                          newdata=data.frame(h_apical = h_apical$x, 
                                                             log_herb_avg = herb_ex[1,1]), 
                                          type="Bertha"),
                    Herbivory = "min")
  avg <- data.frame(h_apical = h_apical$x, 
                    prob.flower = predict(flower.fit, 
                                          newdata=data.frame(h_apical = h_apical$x,
                                                             log_herb_avg = herb_ex[2,1]),
                                          type="Bertha"),
                    Herbivory = "avg")
  max <- data.frame(h_apical = h_apical$x, 
                    prob.flower = predict(flower.fit, 
                                          newdata=data.frame(h_apical = h_apical$x, 
                                                             log_herb_avg = herb_ex[3,1]),
                                          type="Bertha"),
                    Herbivory = "max")
  
  mycurves <- rbind(min, avg, max) 
  mycurves$Herbivory <- as.factor(mycurves$Herbivory)
  
  imgb <- ggplot(mycurves, aes(x = h_apical, y = prob.flower, col = Herbivory, linetype = Herbivory)) +
    geom_line(size = 1.0) + 
    scale_x_continuous(limits = c(0, h_apical$max), expand = c(0, 0)) + 
    scale_color_manual(values = brewer.pal(5, "YlOrRd")[2:4]) +
    xlab("Apical Height (cm)") +
    ylab("Flowering Probability")     
  
  imgb <- imgb + theme(legend.position = "none")
  
  gt <- ggplotGrob(imgt)
  gb <- ggplotGrob(imgb)
  gb <- gtable_add_cols(gb, unit(1, "mm"))
  g <- rbind(gt, gb, size="first")
  g$widths <- unit.pmax(gt$widths, gb$widths)
  
  img <- grid.arrange(g)
  
  detach(obj$pars)
  detach(obj$vars$h_apical)
  
  return(img)
}

renderBudlingDistFit.mwIPM <- function(obj) {
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$vars$h_apical, warn.conflicts = FALSE)

  y <- budling.fit[['BLD1']]$predict(b)
  plotdata <- data.frame(site = "BLD1",
                         x = x,
                         y = y)
  
  y <- budling.fit[['BLD2']]$predict(b)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "BLD2",
                                   x = x,
                                   y = y))
  
  y <- budling.fit[['PWR']]$predict(b)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "PWR",
                                   x = x,
                                   y = y))
  
  y <- budling.fit[['SKY']]$predict(b)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "SKY",
                                   x = x,
                                   y = y))
  
  y <- budling.fit[['YTB']]$predict(b)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "YTB",
                                   x = x,
                                   y = y))
  
  y <- budling.fit[['Bertha']]$predict(b)
  plotdata <- bind_rows(plotdata, 
                        data.frame(site = "Combined",
                                   x = x,
                                   y = y))
  
  plotdata$site <- factor(plotdata$site, levels = c("BLD1", "BLD2", "PWR", "SKY", "YTB", "Combined"))
  
  pb <- plotdata %>% ggplot(aes(x = x, y = y, fill = site, colour = site)) + 
    geom_line() + 
    geom_area(position = "identity", alpha = 0.3)
  
  pb <- pb + ggtitle("Budling Distributions") +
    theme_bw() +
    xlab("Apical Height (cm)") +
    ylab("Probability Density") +
    scale_x_continuous(limits = c(0, 160)) +
    scale_fill_manual(values = c(hue_pal()(5), NA)) + 
    scale_color_manual(values = c(hue_pal()(5), "black")) + 
    labs(colour="Sites", fill="Sites", linetype="Sites") +
    theme(legend.background = element_rect(fill="lightgrey",
                                           size=0.1,
                                           linetype="solid"),
          legend.key.size =  unit(0.2, "in"),
          legend.position = c(0.875, 0.65))
  
  detach(obj$pars)
  detach(obj$vars$h_apical)
  
  return(pb)
}

renderHerbivoryDistFit.mwIPM <- function(obj) {
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$vars$log_herb_avg, warn.conflicts = FALSE)
  
  y <- munched.fit[['BLD1']]$predict(b, justmunch=TRUE)
  plotdata <- data.frame(site = "BLD1",
                         x = x,
                         y = y)
  
  y <- munched.fit[['BLD2']]$predict(b, justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "BLD2",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['PWR']]$predict(b, justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "PWR",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['SKY']]$predict(b, justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "SKY",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['YTB']]$predict(b, justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "YTB",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['Bertha']]$predict(b, justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "Combined",
                                   x = x,
                                   y = y))
  
  plotdata$site <- factor(plotdata$site, levels = c("BLD1", "BLD2", "PWR", "SKY", "YTB", "Combined"))

  pa <- plotdata %>% ggplot(aes(x = x, y = y, fill = site, colour = site)) + 
    geom_line() + 
    geom_area(position = "identity", alpha = 0.3)
  
  pa <- pa + ggtitle("Herbivory Distributions") + 
    theme_bw() +
    xlab("ln(Herbivory Score)") +
    ylab("Probability Density") +
    scale_x_continuous(limits = c(log(0.1), log(6.1))) + 
    scale_fill_manual(values = c(hue_pal()(5), NA)) + 
    scale_color_manual(values = c(hue_pal()(5), "black")) +
    labs(colour="Sites", fill="Sites") +
    theme(legend.background = element_rect(fill="lightgrey",
                                           size=0.1,
                                           linetype="solid"),
          legend.key.size =  unit(0.18, "in"),
          legend.position = c(0.885, 0.72))
  
  detach(obj$pars)
  detach(obj$vars$log_herb_avg)
  
  return(pa)
}