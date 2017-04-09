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
setHerbivoryMatrix <- function(obj, dist.herb, update, perturb) UseMethod("setHerbivoryMatrix")
setSeedlingRecruitmentMatrix <- function(obj, update, perturb) UseMethod("setSeedlingRecruitmentMatrix")

# Compute kernels from matrices
computeSexualKernel <- function(obj, update, perturb) UseMethod("computeSexualKernel")
computeClonalKernel <- function(obj, update, perturb) UseMethod("computeClonalKernel")
computeFullKernel <- function(obj) UseMethod("computeFullKernel")

# Setup or compute MPM/IPM
computeMPM <- function(obj) UseMethod("computeMPM")
setSite <- function(obj, site, compute) UseMethod("setSite")
bootIPM <- function(obj) UseMethod("bootIPM")

# Analyze growth rate
analyzeGrowthRate <- function(obj) UseMethod("analyzeGrowthRate")
analyzeStandard <- function(obj) UseMethod("analyzeStandard")
analyzeParameters <- function(obj, compute, saveresults, params) UseMethod("analyzeParameters")

# Renderers
renderFloweringFit <- function(obj) UseMethod("renderFloweringFit") 
renderBudlingDistFit <- function(obj) UseMethod("renderBudlingDistFit")
renderHerbivoryDistFit <- function(obj) UseMethod("renderHerbivoryDistFit")

# Helpers & Globals

glmerCtrl <- glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun=50000))

# Keeping these commented here as a record
# tfunc <- function(x, y)
# perturbTrans <- function(pars, perturb)

# Constructor ----------------------

mwIPM <- function(x = list()) {
  #
  # Specify min and max of variables in list form.
  # Note:  Had to move this AFTER bootIPM or it would'nt load - ?????
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
    saveresults = (length(list.files(mwROOT("data","calculated"), pattern=".RData")) == 0)
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
                      MPM = NA,
                      analysis = list(standard = NA,
                                      parameters = NA)))
  
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
      obj <- setSeedsPerPodConst(obj)
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
                                        !is.na(log_herb_avg),
                                        !is.na(fec.flower))
    
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
    h_apical <- (obj$data %>% filter(seedling == 1, !is.na(h_apical)))$h_apical
    
    cat("Computing seedling distribution fit...")
    f1 <- fitdist(h_apical, "lnorm")
    cat("done!\n")
    
    seedling.fit <- vector("list", 2)
    seedling.fit[[1]] <- f1
    
    # NOTE:  Doesn't need to be eval(parse(text = ...))!!!! Just make the function, silly...  will change later.
    seedling.fit[[2]] <- eval(parse(
      text = 
        "function(x, pars, perturb = rep(0,2)) {
           pars <- perturbTrans(pars, perturb)
           N <- length(x)
           dx <- x[2]-x[1]
           y <- rep(0, N-1)
           for (j in 1:(N-1)) {
             y[j] = plnorm(x[j+1], pars[1], pars[2]) - plnorm(x[j], pars[1], pars[2])
           }
           y <- y/(dx*sum(y))
         }"
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
    mdls <- c("norm", "norm", "norm", "norm", "weibull", "weibull")
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
            "function(x, pars, perturb = rep(0,2)) {
               pars <- pars + perturb
               N <- length(x)
               dx <- x[2]-x[1]
               y <- rep(0, N-1)
               for (j in 1:(N-1)) {
                 y[j] = p%s(x[j+1], pars[1], pars[2]) - p%s(x[j], pars[1], pars[2])
               }
               y <- y/(dx*sum(y))
            }",
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$distname
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
    mdls <- c("lnorm", "lnorm", "lnorm", "lnorm", "lnorm", "gamma")
    munched.fit <- vector('list', 6)
    names(munched.fit) <- sites
    
    for (i in 1:6) {
      if (i == 1) { # Bertha
        thissite <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(munched))
        functext <- "function(x, pars, perturb = rep(0,3), justmunch = FALSE) {
                       pars[1] <- pars[1] + perturb[1]\n
                       pars[2:3] <- perturbTrans(pars[2:3], perturb[2:3])\n"
      } else { # Sites
        thissite <- obj$data %>% filter(site == sites[i],
                                        !is.na(h_apical),
                                        !is.na(munched))
        functext <- "function(x, pars, perturb = rep(0,3), justmunch = FALSE) {
                       pars <- pars + perturb\n"
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
        eval(parse(text = paste0(functext,
          sprintf("N <- length(x)
                   dx <- x[2]-x[1]
                   z <- exp(x)-0.1
                   y <- rep(0, N)
                   for (j in 1:(N-1)) {
                     y[j] = p%s(z[j+1], pars[2], pars[3]) - p%s(z[j], pars[2], pars[3])
                   }
                   y <- (pars[1]/(dx*sum(y)))*y
                   if (!justmunch) {
                     y[1] <- y[1] + (1-pars[1])/dx
                   }
                   y <- y[1:(N-1)]
                 }",
                 munched.fit[[i]][[1]]$distname,
                 munched.fit[[i]][[1]]$distname
                 )
        )))
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

    budlings.per.stem.fit <- vector("list", 3)
    budlings.per.stem.fit[[1]] <- exp.model
    budlings.per.stem.fit[[2]] <- merged$site
    budlings.per.stem.fit[[3]] <- function (x, pars, perturb = rep(0,2)) {
      pars <- pars + perturb
      y <- exp(pars[1] + pars[2]*x)
    }
    names(budlings.per.stem.fit) <- c("fit","site","predict")

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

  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }
  
  return(obj)
}

setSurvivalMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Survival
  # N x N^2
  obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]
  
  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }
  
  return(obj)
}

setGrowthMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,5)) {
  # Growth
  # N x N^2

  N <- obj$N
  Mu <- c(outer(obj$vars$h_apical$x, obj$vars$log_herb_avg$x, function(x,y) {predict(obj$pars$growth.fit, newdata = data.frame(h_apical = x, log_herb_avg = y), type=obj$site, perturb=perturb[1:4])}))
  G <- matrix(rep(0, N^3), nrow=N)
  for (i in (1:N^2)) {
    for (j in 1:N) {
      G[j,i] = pnorm(obj$vars$h_apical.next$b[j+1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]]+perturb[5]) - pnorm(obj$vars$h_apical.next$b[j], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]]+perturb[5])
    }
    G[1,i] = G[1,i] + pnorm(obj$vars$h_apical.next$b[1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]]+perturb[5])
    G[N,i] = G[N,i] + pnorm(obj$vars$h_apical.next$b[N+1], Mu[i], obj$pars$growth.fit$pars$sd[[obj$site]]+perturb[5], lower.tail=FALSE)
    G[,i] <- G[,i]/obj$vars$h_apical.next$dx
  }
  obj$matrices$G <- G
  
  # Check for eviction
  # colSums(G%*%H*dx.h_apical.next*dx.log_herb_avg)
  
  # Plot growth
  # image.plot(h_apical, h_apical.next, t(G%*%H), col=topo.colors(100))
  
  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }
  
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
  
  if (update) {
    obj %<>% computeSexualKernel()
    obj <- computeMPM(obj)
  }
  
  return(obj)
}

setHerbivoryMatrix.mwIPM <- function(obj, dist.herb = NA, update = TRUE, perturb = rep(0,3)) {
  # Herbivory matrix
  # N x N^2
  
  N <- obj$N
  if (any(is.na(dist.herb))) {
    dist.herb <- obj$pars$munched.fit[[obj$site]]$predict(obj$vars$log_herb_avg$b,
                                                          c(obj$pars$munched.fit[[obj$site]]$pmunch, 
                                                            obj$pars$munched.fit[[obj$site]]$fit$estimate),
                                                          perturb)
  }

  H = matrix(rep(0,N^3), nrow = N^2)
  for (i in 1:N) {
    for (I in 1:N) {
      H[I+N*(i-1), I] = dist.herb[i]
    }
  }

  obj$pars$dist.herb <- dist.herb
  obj$matrices$H <- H
  
  if (update) {
    obj %<>% computeSexualKernel(update = FALSE) %>%
             computeClonalKernel(update = FALSE) %>%
             computeFullKernel()
  }
  
  return(obj)
}

setSeedlingRecruitmentMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,2)) {
  # N x N
  
  R <- t(t(obj$pars$seedling.fit$predict(obj$vars$h_apical$b,
                                         obj$pars$seedling.fit$fit$estimate,
                                         perturb)))%*%t(rep(1,obj$N))
  
  obj$matrices$R <- R
  
  if (update) {
    obj %<>% computeSexualKernel()
    obj <- computeMPM(obj)
  }
  
  return(obj)
}

# Kernels ------------------------------

computeSexualKernel.mwIPM <- function(obj, update = TRUE, perturb = rep(0,2)) {
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)
  
  Ks <- (seedling.emergence[[obj$site]]+perturb[1])*(seeds.per.pod+perturb[2])*R%*%(P*G*S*F)%*%H*log_herb_avg$dx*h_apical.next$dx*h_apical$dx
  # image.plot(h_apical, h_apical, t(Ks), col=topo.colors(100))
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$kernels$Ks <- Ks
  
  if (update) {
    obj %<>% computeFullKernel()
  }
  
  return(obj)
}

computeClonalKernel.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  ## Budling Recruitment (same as Kc)
  # perturb = c(budlings.per.stem, budlings)
  
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)
  
  mean.log_herb_avg <- sum(dist.herb*log_herb_avg$x)*log_herb_avg$dx
  # mean.buds.per.stem <- exp(predict(budlings.per.stem.fit, new = data.frame(log_herb_mean = mean.log_herb_avg)))
  
  mean.buds.per.stem <- budlings.per.stem.fit$predict(mean.log_herb_avg,
                                                      budlings.per.stem.fit$fit$coefficients,
                                                      perturb[1:2])
  Kc <- t(t(mean.buds.per.stem*budling.fit[[obj$site]]$predict(h_apical$b,
                                                               budling.fit[[obj$site]]$fit$estimate,
                                                               perturb[3:4])))%*%t(rep(1,obj$N))*h_apical$dx
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$kernels$Kc <- Kc
  
  if (update) {
    obj %<>% computeFullKernel()
  }
  
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
      
      obj %<>% setVars()
      
      obj %<>% 
        setFloweringMatrix(update = FALSE) %>% 
        setSurvivalMatrix(update = FALSE) %>% 
        setGrowthMatrix(update = FALSE) %>% 
        setPodsMatrix(update = FALSE) %>%
        setHerbivoryMatrix(update=FALSE) %>% 
        setSeedlingRecruitmentMatrix(update = FALSE)
      
      obj %<>%
        computeSexualKernel(update = FALSE) %>% 
        computeClonalKernel(update = FALSE) %>%
        computeFullKernel()
      
      obj %<>% computeMPM()
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
  # mean.buds.per.stem <- exp(predict(budlings.per.stem.fit, new = data.frame(log_herb_mean = mean.log_herb_avg)))
  mean.buds.per.stem <- budlings.per.stem.fit$predict(mean.log_herb_avg,
                                                      budlings.per.stem.fit$fit$coefficients)
  
  Kss <- seedling.emergence[[obj$site]]*seeds.per.pod*(P*G*S*F)%*%H*log_herb_avg$dx*h_apical.next$dx*h_apical$dx
  alpha <- sum(Kss%*%t(t(seedling.fit$predict(h_apical$b, seedling.fit$fit$estimate))))
  beta <- sum(Kss%*%t(t(budling.fit[[obj$site]]$predict(h_apical$b,
                                                        budling.fit[[obj$site]]$fit$estimate))))
  pems <- mean.buds.per.stem
  pemb <- mean.buds.per.stem
  
  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)
  
  obj$MPM <- matrix(c(alpha, beta, pems, pemb), nrow = 2, byrow=TRUE)
  
  return(obj)
}

analyzeGrowthRate.mwIPM <- function(obj) {
  return(Re(eigen(obj$kernels$K)$values[1]))
}

analyzeStandard.mwIPM <- function(obj) {
  K <- obj$kernels$K
  lam <- analyzeGrowthRate(obj)
  w.eigen <- Re(eigen(K)$vectors[,1])
  stable.dist <- w.eigen/(obj$vars$h_apical$dx*sum(w.eigen))
  v.eigen <- Re(eigen(t(K))$vectors[,1])
  repro.val <- v.eigen/v.eigen[1]
  v.dot.w <- sum(stable.dist*repro.val)*obj$vars$h_apical$dx
  sens <- outer(repro.val,stable.dist)/v.dot.w
  elas <- matrix(as.vector(sens)*as.vector(K)/lam,nrow=ipm$N)
  obj$analysis$standard <- list(lambda = lam,
                                stable.dist = stable.dist,
                                repro.val = repro.val,
                                sens = sens,
                                elas = elas)
  return(obj)
}

analyzeParameters.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, params = "all") {
  if (!file.exists(mwROOT("data","calculated","parameterAnalysis.RData")) | (compute)) {
    # Initialize to empty data frame
    analysis <- tbl_df(data.frame(sensitivity = NULL,
                                  pars = NULL,
                                  type = NULL,
                                  name = NULL))
    
    if (("all" %in% params) | ("Flowering" %in% params)) {
      cat("Analyzing flowering fit...")
      flowering_func <- function(x) {obj %>% setFloweringMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(flowering_func,rep(0,4)),
                          pars = obj$pars$flower.fit$pars$unscaled["Bertha",],
                          type = as.character("Flowering"),
                          name = c("(Intercept)", 
                                   "h_apical", 
                                   "log_herb_avg", 
                                   "h_apical:log_herb_avg")
        )
        )
      )
      cat("done!\n")
    }
    
    if (("all" %in% params) | ("Survival" %in% params)) {
      cat("Analyzing survival fit...")
      survival_func <- function(x) {obj %>% setSurvivalMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(survival_func, rep(0,4)),
                          pars = obj$pars$surv.fit$pars$unscaled["Bertha",],
                          type = as.character("Survival"),
                          name = c("(Intercept)", 
                                   "h_apical", 
                                   "log_herb_avg", 
                                   "h_apical:log_herb_avg")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Growth" %in% params)) {
      cat("Analyzing growth fit...")
      growth_func <- function(x) {obj %>% setGrowthMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(growth_func, rep(0,5)),
                          pars = c(obj$pars$growth.fit$pars$unscaled["Bertha",],
                                   obj$pars$growth.fit$pars$sd["Bertha"]),
                          type = as.character("Growth"),
                          name = c("(Intercept)", 
                                   "h_apical", 
                                   "log_herb_avg", 
                                   "h_apical:log_herb_avg", 
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Pods" %in% params)) {    
      cat("Analyzing pods fit...")
      pods_func <- function(x) {obj %>% setPodsMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(pods_func, rep(0,4)),
                          pars = obj$pars$pods.fit$pars$unscaled["Bertha",],
                          type = as.character("Pods"),
                          name = c("(Intercept)", 
                                   "h_apical", 
                                   "log_herb_avg", 
                                   "h_apical:log_herb_avg")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Seedlings" %in% params)) {    
      cat("Analyzing seedling distribution...")
      seedling_func <- function(x) {obj %>% setSeedlingRecruitmentMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(seedling_func, rep(0,2)),
                          pars = tfunc(x = obj$pars$seedling.fit$fit$estimate, y = c(0,0)),
                          type = as.character("Seedlings"),
                          name = c("mean", 
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }
    
    if (("all" %in% params) | ("Sexual" %in% params)) {
      cat("Analyzing seedling emergence and seeds-per-pod...")
      sexual_func <- function(x) {obj %>% computeSexualKernel(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(sexual_func, rep(0,2)),
                          pars = c(obj$pars$seedling.emergence["Bertha"], 
                                   obj$pars$seeds.per.pod),
                          type = as.character("Sexual"),
                          name = c("seedling.emergence", 
                                   "seeds.per.pod")
        )
        )
      )
      cat("done!\n")
    }
    
    if (("all" %in% params) | ("Clonal" %in% params)) {
      cat("Analyzing budlings-per-stem and budling distribution...")
      clonal_func <- function(x) {obj %>% computeClonalKernel(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(clonal_func, rep(0,4)),
                          pars = c(obj$pars$budlings.per.stem.fit$fit$coefficients,
                                   obj$pars$budling.fit$Bertha$fit$estimate),
                          type = c("Clonal", "Clonal", "Budlings", "Budlings"),
                          name = c("(Intercept)", 
                                   "log_herb_avg",
                                   "mean",
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Herbivory" %in% params)) {    
      cat("Analyzing herbivory distribution...")
      herbivory_func <- function(x) {obj %>% setHerbivoryMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = grad(herbivory_func, rep(0,3)),
                          pars = c(obj$pars$munched.fit$Bertha$pmunch,
                                   tfunc(x = obj$pars$munched.fit$Bertha$fit$estimate, y = c(0,0))),
                          type = c("Herbivory"),
                          name = c("pmunch", 
                                   "mean",
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }
    
    analysis %<>% filter(pars != 0) %>%
      mutate(lambda = ipm %>% analyzeGrowthRate(),
             elasticity = sensitivity/(lambda/pars))
    analysis$type <- factor(analysis$type)
    
    if (saveresults) {
      save(analysis, file = mwROOT("data","calculated","parameterAnalysis.RData"))
    }
  } else {
    load(mwROOT("data","calculated","parameterAnalysis.RData"))
  }
  
  obj$analysis$parameters <- analysis
  
  return(obj)
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
  
  y <- munched.fit[['BLD1']]$predict(b, 
                                     c(munched.fit[['BLD1']]$pmunch, 
                                       munched.fit[['BLD1']]$fit$estimate),
                                     justmunch=TRUE)
  plotdata <- data.frame(site = "BLD1",
                         x = x,
                         y = y)
  
  y <- munched.fit[['BLD2']]$predict(b,
                                     c(munched.fit[['BLD2']]$pmunch, 
                                       munched.fit[['BLD2']]$fit$estimate),
                                     justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "BLD2",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['PWR']]$predict(b,
                                    c(munched.fit[['PWR']]$pmunch, 
                                      munched.fit[['PWR']]$fit$estimate),
                                    justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "PWR",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['SKY']]$predict(b,
                                    c(munched.fit[['SKY']]$pmunch, 
                                      munched.fit[['SKY']]$fit$estimate),
                                    justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "SKY",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['YTB']]$predict(b,
                                    c(munched.fit[['YTB']]$pmunch, 
                                      munched.fit[['YTB']]$fit$estimate),
                                    justmunch=TRUE)
  plotdata <- bind_rows(plotdata,
                        data.frame(site = "YTB",
                                   x = x,
                                   y = y))
  
  y <- munched.fit[['Bertha']]$predict(b,
                                       c(munched.fit[['Bertha']]$pmunch, 
                                         munched.fit[['Bertha']]$fit$estimate),
                                       justmunch=TRUE)
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

# Helpers  ------------------------------

tfunc <- function(x, y) {
  # In all cases that we care about (i.e. Bertha), log-normal is the only one requiring
  #   a transformation of the perturbation.
  # In the log-normal case, x is (location, scale), and y is (mean, sd).
  c(F1 = exp(x[1] + 0.5*x[2]*x[2]) - y[1],
    F2 = sqrt(exp(2*x[1] + x[2]*x[2])*(exp(x[2]*x[2])-1)) - y[2])
}

perturbTrans <- function(pars, perturb = rep(0,2)) {
  # Let's first get the mean and sd from the location and scale.
  MSD <- tfunc(x = pars, y = c(0,0))
  
  # Now we perturb the mean and sd and use multiroot to find the new
  #   location and scale, using the old location and scale as the initial state.
  tpars <- multiroot(tfunc, pars, y = MSD+perturb)$root
}