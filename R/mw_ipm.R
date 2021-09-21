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

#' Create herbivory matrix.
#'
#' @param obj A mwIPM model object.
#' @param dist.herb Herbivory distribution.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#' @param distpars Determines in which scale (log or linear) perturbation is occurring.
#'
#' @return A mwIPM model object.
#' @export
#'
#' @importFrom magrittr %<>%
setHerbivoryMatrix <- function(obj, dist.herb, update, perturb, distpars) UseMethod("setHerbivoryMatrix")

setSeedlingRecruitmentMatrix <- function(obj, update, perturb) UseMethod("setSeedlingRecruitmentMatrix")

# Compute kernels from matrices
computeSexualKernel <- function(obj, update, perturb) UseMethod("computeSexualKernel")
computeClonalKernel <- function(obj, update, perturb) UseMethod("computeClonalKernel")
computeFullKernel <- function(obj) UseMethod("computeFullKernel")

# Setup or compute MPM/IPM
computeMPM <- function(obj) UseMethod("computeMPM")

#' Set site and recompute kernel.
#'
#' @param obj A mwIPM model object.
#' @param site Site ("Bertha", "BLD1", "BLD2", "GET", "GRN", "PWR", "SKY", "YTB")
#' @param compute Recompute or load from cache?
#'
#' @return A mwIPM model object.
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @examples
#' ipm %<>% setSite("BLD1")
setSite <- function(obj, site, compute) UseMethod("setSite")

#' Bootstraps over the stems.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object.
#' @export
#' @import dplyr
#' @importFrom magrittr %>% %<>%
#'
#' @examples
#' ipm %<>% bootIPM()
bootIPM <- function(obj) UseMethod("bootIPM")

# Analyses

#' Compute the population growth rate.
#'
#' @param obj A mwIPM model object.
#'
#' @return The population growth rate
#' @export
#'
#' @examples
#' lambda <- ipm %>% analyzeGrowthRate()
analyzeGrowthRate <- function(obj) UseMethod("analyzeGrowthRate")

#' Compute standard IPM population level metrics, such as the population
#' growth rate, left and right eigenvectors, and sensitivity and elasticit
#' matrices.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object, with the results of the computation stored
#' in a list.
#' @export
#'
#' @examples
#' ipm %<>% analyzeStandard()
#' str(ipm$analysis$standard)
analyzeStandard <- function(obj) UseMethod("analyzeStandard")

#' Sensitivity analysis on parameters.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param params Which parameter to perturb?  One of "all", "Flowering", ...
#' @param distpars Determines in which scale (log or linear) perturbation is occurring.
#'
#' @return A mwIPM model object, with the results of the computation stored
#' in a list.
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#'
#' @examples
#' ipm %<>% analyzeParameters()
#' str(ipm$analysis$parameters)
analyzeParameters <- function(obj, compute, saveresults, params, distpars) UseMethod("analyzeParameters")

# Renderers
renderFloweringFit <- function(obj) UseMethod("renderFloweringFit")
renderBudlingDistFit <- function(obj) UseMethod("renderBudlingDistFit")
renderHerbivoryDistFit <- function(obj) UseMethod("renderHerbivoryDistFit")

# Helpers & Globals

glmerCtrl <- lme4::glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun=50000))
mwCache <- file.path(rappdirs::app_dir("Milkweed", "LaMar")$data(), "calculated")
options(milkweed.cache = mwCache) # Store this for outside use in vignettes

# Keeping these commented here as a record
# ParsToMoms <- function(x)
# MomsToPars <- function(y)
# perturbTrans <- function(pars, perturb)

# Constructor ----------------------

#' Construcor method for mwIPM S3 class.
#'
#' @param x A list of options.
#'
#' @return A mwIPM model object.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @examples
#' ipm <- mwIPM()
mwIPM <- function(x = list()) {
  x$all_sites = c("Bertha", levels(stemdata$site))
  if (all(names(x) != "N")) {
    x$N = 50
  }
  if (all(names(x) != "data")) {
    x$data <- stemdata %>% filter(site %in% x$all_sites)
  }
  if (all(names(x) != "site")) {
    site = "Bertha"
    x$site <- NA
  } else if (!(x$site %in% x$all_sites)) {
    stop(paste(x$site, "is not a valid site!  Please choose one of ", paste0(toString(x$all_sites), ".")))
  }
  if (all(names(x) != "compute")) {
    compute = FALSE
  } else if (!is.logical(x$compute)) {
    stop(paste("Variable 'compute' must be TRUE or FALSE."))
  } else {
    compute <- x$compute
  }
  if (all(names(x) != "saveresults")) {
    saveresults = (length(list.files(mwCache, pattern=".RData")) == 0)
  } else if (!is.logical(x$saveresults)) {
    stop(paste("Variable 'saveresults' must be TRUE or FALSE."))
  } else {
    saveresults <- x$saveresults
  }
  if (all(names(x) != "mdlargs")) {
    x$mdlargs <- list(method = 'linear',
                      input = 'full')
  } else {
    if (all(names(x$mdlargs) != "method")) {
      x$mdlargs$method = "linear"
    }
    if (all(names(x$mdlargs) != "input")) {
      x$mdlargs$input = "full"
    }
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
                                  herb_avg = list(min = NA,
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

#' Bootstraps over the stems.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object.
#' @export
#' @import dplyr
#' @importFrom magrittr %>% %<>%
#'
#' @examples
#' ipm %<>% bootIPM()
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

#' Sets model variables.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
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

  # herb_avg
  herb_avg <- list(min = 0.0,
                   max = 6.0)
  herb_avg$b = herb_avg$min + c(0:N)*(herb_avg$max - herb_avg$min)/N # boundary points
  herb_avg$x = 0.5*(herb_avg$b[1:N] + herb_avg$b[2:(N+1)])
  herb_avg$dx = herb_avg$b[2] - herb_avg$b[1] # class size

  obj$vars = list(h_apical = h_apical,
                  h_apical.next = h_apical.next,
                  herb_avg = herb_avg)

  return(obj)
}

#' Sets model parameters.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
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

#' Sets the number of seeds per pod parameter.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
setSeedsPerPodConst.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"seedsPerPodConst.RData")) | (compute)) {
    seeds_per_pod_data <- seeddata

    cat("Computing seeds per pod constant...")
    seeds.per.pod = mean(seeds_per_pod_data$total_seed)
    cat("done!\n")

    if (saveresults) {
      save(seeds.per.pod, file = file.path(mwCache,"seedsPerPodConst.RData"))
    }
  } else {
    load(file.path(mwCache,"seedsPerPodConst.RData"))
  }

  obj$pars$seeds.per.pod <- seeds.per.pod

  return(obj)
}

#' Sets the seedling emergence parameter.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setSeedlingEmergenceConst.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"seedlingEmergenceConst.RData")) | (compute)) {
    seedling.emergence <- rep(NA, length(obj$all_sites))
    names(seedling.emergence) <- obj$all_sites

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

    # Sites - really need to map this one
    data_gp <- obj$data %>% group_by(year, site) %>%
                            summarize(N_seedlings = sum(seedling, na.rm=T),
                                      N_seeds = seeds.per.pod*sum(N_pods, na.rm=T))

    data13_14 <- data_gp %>% filter(year %in% 2013:2014) %>%
                             group_by(site) %>%
                             summarize(emergence = last(N_seedlings)/first(N_seeds))

    data14_15 <- data_gp %>% filter(year %in% 2014:2015) %>%
                             group_by(site) %>%
                             summarize(emergence = last(N_seedlings)/first(N_seeds))

    data15_16 <- data_gp %>% filter(year %in% 2015:2016) %>%
      group_by(site) %>%
      summarize(emergence = last(N_seedlings)/first(N_seeds))

    data16_17 <- data_gp %>% filter(year %in% 2016:2017) %>%
      group_by(site) %>%
      summarize(emergence = last(N_seedlings)/first(N_seeds))

    fulldat <- bind_rows(data13_14, data14_15, data15_16, data16_17)

    seedling.emergence[2:length(obj$all_sites)] <- (fulldat %>% group_by(site) %>% summarize(emergence = mean(emergence, na.rm = TRUE)))$emergence
    # GRN is at 0 - set to Bertha
    seedling.emergence["GRN"] <- seedling.emergence["Bertha"]
    cat("done!\n")

    if (saveresults) {
      save(seedling.emergence, file = file.path(mwCache,"seedlingEmergenceConst.RData"))
    }
  } else {
    load(file.path(mwCache,"seedlingEmergenceConst.RData"))
  }

  obj$pars$seedling.emergence <- seedling.emergence

  return(obj)
}

# Fits ------------------------------

## Vital rates

#' Sets the flowering mixed-model.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setFloweringFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"flowerFit.RData")) | (compute)) {
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(herb_avg),
                                        !is.na(fec.flower))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, herb_avg),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing flowering fit...")
    flower.mdl <- lme4::glmer(fec.flower ~ h_apical*herb_avg + (1|site/transect)+(h_apical+herb_avg|year),
                              data=metadata_sc,
                              nAGQ=1,
                              family=binomial(),
                              control=glmerCtrl)
    cat("done!\n")

    flower.fit <- mwMod(list(mdl = flower.mdl,
                             vars = c("h_apical", "herb_avg"),
                             scaled = list(h_apical = scale(metadata_usc$h_apical),
                                           herb_avg = scale(metadata_usc$herb_avg))))

    # Check parameters
    cat("Checking parameters:\n")
    checkPars(flower.fit)

    if (saveresults) {
      save(flower.fit, file = file.path(mwCache,"flowerFit.RData"))
    }
  } else {
    load(file.path(mwCache,"flowerFit.RData"))
  }

  obj$pars$flower.fit <- flower.fit

  return(obj)
}

#' Sets the survival mixed-model.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setSurvivalFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"survivalFit.RData")) | (compute)) {
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        !is.na(surv))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, herb_avg),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing survival fit...")
    surv.mdl <- lme4::glmer(surv ~ h_apical + herb_avg + (h_apical|site/transect) + (h_apical|year),
                            data=metadata_sc,
                            family=binomial(),
                            nAGQ=1,
                            control=glmerCtrl)
    cat("done!\n")

    surv.fit <- mwMod(list(mdl = surv.mdl,
                           vars = c("h_apical", "herb_avg"),
                           scaled = list(h_apical = scale(metadata_usc$h_apical),
                                         herb_avg = scale(metadata_usc$herb_avg))))

    # Check parameters
    cat("Checking parameters:\n")
    checkPars(surv.fit)

    if (saveresults) {
      save(surv.fit, file = file.path(mwCache,"survivalFit.RData"))
    }
  } else {
    load(file.path(mwCache,"survivalFit.RData"))
  }

  obj$pars$surv.fit <- surv.fit

  return(obj)
}

#' Sets the growth mixed-model.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setGrowthFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"growthFit.RData")) | (compute)) {
    # We only want stems that flowered and survived
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        surv == 1)

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, h_apical.next, herb_avg),
                                              .funs = funs(as.numeric(scale(.))))


    cat("Computing growth fit...")
     growth.mdl <- lme4::lmer(h_apical.next ~ h_apical*herb_avg + (1|site/transect) + (h_apical|year),
                              data=metadata_sc,
                              REML=T)
    cat("done!\n")

    growth.fit <- mwMod(list(mdl = growth.mdl,
                             vars = c("h_apical", "herb_avg"),
                             scaled = list(h_apical = scale(metadata_usc$h_apical),
                                           h_apical.next = scale(metadata_usc$h_apical.next),
                                           herb_avg = scale(metadata_usc$herb_avg))))
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(growth.fit)

    if (saveresults) {
      save(growth.fit, file = file.path(mwCache,"growthFit.RData"))
    }
  } else {
    load(file.path(mwCache,"growthFit.RData"))
  }

  obj$pars$growth.fit <- growth.fit

  return(obj)
}

#' Sets the pods mixed-model.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setPodsFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"podsFit.RData")) | (compute)) {
    # We only want stems that flowered, survived and have data for pods
    metadata_usc <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        surv == 1,
                                        !is.na(N_pods))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical.next, herb_avg),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing pods fit...")
    pods.mdl <- lme4::glmer(N_pods ~ h_apical.next + (h_apical.next|site/transect) + (herb_avg|year),
                            data=metadata_sc,
                            nAGQ=1,
                            family=poisson(link = "log"),
                            control=glmerCtrl)
    cat("done!\n")

    pods.fit <- mwMod(list(mdl = pods.mdl,
                           vars = c("h_apical.next", "herb_avg"),
                           scaled = list(h_apical.next = scale(metadata_usc$h_apical.next),
                                         herb_avg = scale(metadata_usc$herb_avg))))
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(pods.fit)

    if (saveresults) {
      save(pods.fit, file = file.path(mwCache,"podsFit.RData"))
    }
  } else {
    load(file.path(mwCache,"podsFit.RData"))
  }

  obj$pars$pods.fit <- pods.fit

  return(obj)
}

## Distributions

#' Sets the seedling distribution (normal distribution).
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setSeedlingDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"seedlingDistFit.RData")) | (compute)) {
    h_apical <- (obj$data %>% filter(seedling == 1, !is.na(h_apical)))$h_apical

    cat("Computing seedling distribution fit...")
    f1 <- fitdistrplus::fitdist(h_apical, "norm")
    cat("done!\n")

    seedling.fit <- vector("list", 2)
    seedling.fit[[1]] <- f1

    seedling.fit[[2]] <- function(x, pars, perturb = rep(0,2)) {
           pars <- perturbTrans(pars, perturb)
           N <- length(x)
           dx <- x[2]-x[1]
           y <- rep(0, N-1)
           for (j in 1:(N-1)) {
             y[j] = pnorm(x[j+1], pars[1], pars[2]) - pnorm(x[j], pars[1], pars[2])
           }
           y[1] <- y[1] + pnorm(x[1], pars[1], pars[2])
           y[N-1] <- y[N-1] + pnorm(x[N], pars[1], pars[2], lower.tail = FALSE)
           y <- y/dx
         }
    names(seedling.fit) <- c("fit", "predict")

    if (saveresults) {
      save(seedling.fit, file = file.path(mwCache,"seedlingDistFit.RData"))
    }
  } else {
    load(file.path(mwCache,"seedlingDistFit.RData"))
  }

  obj$pars$seedling.fit <- seedling.fit

  return(obj)
}

#' Sets the budling distribution.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setBudlingDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"budlingDistFit.RData")) | (compute)) {
    sites <- obj$all_sites
    num_sites <- length(obj$all_sites)
    mdls <- rep("norm", num_sites)
    budling.fit <- vector('list', num_sites)
    names(budling.fit) <- sites

    for (i in 1:num_sites) {
      if (i == 1) { # Bertha
        h_apical <- (obj$data %>% filter(seedling == 0,
                                         !is.na(h_apical)))$h_apical
      } else { # Sites
        h_apical <- (obj$data %>% filter(seedling == 0,
                                         site == sites[i],
                                         !is.na(h_apical)))$h_apical
      }

      cat("Computing budling distribution fit for", sites[i], "...")
      f0 <- fitdistrplus::fitdist(h_apical, mdls[i])
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
               y[1] <- y[1] + p%s(x[1], pars[1], pars[2])
               y[N-1] <- y[N-1] + p%s(x[N], pars[1], pars[2], lower.tail = FALSE)
               y <- y/dx
            }",
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$distname,
            budling.fit[[i]][[1]]$distname
          )
        ))
      names(budling.fit[[i]]) <- c("fit", "predict")
    }

    if (saveresults) {
      save(budling.fit, file = file.path(mwCache,"budlingDistFit.RData"))
    }
  } else {
    load(file.path(mwCache,"budlingDistFit.RData"))
  }

  obj$pars$budling.fit <- budling.fit

  return(obj)
}

#' Sets the herbivory distribution.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
setHerbivoryDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"herbivoryDistFit.RData")) | (compute)) {
    num_sites <- length(obj$all_sites)

    mdls <- rep("lnorm", num_sites)
    names(mdls) <- obj$all_sites

    munched.fit <- vector('list', num_sites)
    names(munched.fit) <- obj$all_sites

    for (i in 1:num_sites) {
      if (i == 1) { # Bertha
        thissite <- obj$data %>% filter(!is.na(h_apical),
                                        !is.na(munched))
      } else { # Sites
        thissite <- obj$data %>% filter(site == obj$all_sites[i],
                                        !is.na(h_apical),
                                        !is.na(munched))
      }
      functext <-
        sprintf(
          paste0("function(x, pars, perturb = rep(0,3), distpars = FALSE, justmunch = FALSE) {
               pars[1] <- pars[1] + perturb[1]
               pars[2:3] <- perturbTrans(pars[2:3], perturb[2:3], type = ifelse(distpars, '%s', 'ident'))\n"),
          mdls[i]
        )
      pmunch <- sum(thissite$munched == 1)/nrow(thissite)
      herb_avg <- (thissite %>% filter(munched == 1))$herb_avg

      cat("Computing herbivory distribution fit for", obj$all_sites[i], "...")
      f0 <- fitdistrplus::fitdist(herb_avg, mdls[i])
      cat("done!\n")

      munched.fit[[i]] <- vector("list", 3)
      munched.fit[[i]][[1]] <- f0
      munched.fit[[i]][[2]] <- pmunch
      munched.fit[[i]][[3]] <-
        eval(parse(text = paste0(functext,
          sprintf("N <- length(x)
                   dx <- x[2]-x[1]
                   y <- rep(0, N-1)
                   for (j in 1:(N-1)) {
                     y[j] = p%s(x[j+1], pars[2], pars[3]) - p%s(x[j], pars[2], pars[3])
                   }
                   y[1] <- y[1] + p%s(x[1], pars[2], pars[3])
                   y[N-1] <- y[N-1] + p%s(x[N], pars[2], pars[3], lower.tail = FALSE)
                   y <- pars[1]*y
                   if (!justmunch) {
                     y[1] <- y[1] + (1-pars[1])
                   }
                   y <- y/dx
                 }",
                 munched.fit[[i]][[1]]$distname,
                 munched.fit[[i]][[1]]$distname,
                 munched.fit[[i]][[1]]$distname,
                 munched.fit[[i]][[1]]$distname
                 )
        )))
      names(munched.fit[[i]]) <- c("fit", "pmunch", "predict")
    }

    if (saveresults) {
      save(munched.fit, file = file.path(mwCache,"herbivoryDistFit.RData"))
    }
  } else {
    load(file.path(mwCache,"herbivoryDistFit.RData"))
  }

  obj$pars$munched.fit <- munched.fit

  return(obj)
}

#' Set budlings-per-stem model.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param update Update dependencies?
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
#' @import dplyr
setBudlingsPerStemFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path(mwCache,"budlingsPerStemFit.RData")) | (compute)) {
    # Note that we are using stemdata and not obj$data!!!!  We can refer to obj$data once we have
    #  full integrated and cleaned data for 2017.  This should be the ONLY place where we use
    #  stemdata in place of obj$data.
    data_gp <- stemdata %>%
      group_by(year, transect) %>%
      summarize(N_seedlings = sum(seedling, na.rm=T),
                N_total = sum(aliveJune, na.rm=T),
                N_budlings = N_total - N_seedlings,
                herb_mean = mean(herb_avg, na.rm=T),
                site = first(site)) %>%
      ungroup(year, transect) %>%
      ##NOTE: these 4 bind_rows() calls add in explicit 0s where transects died off, resulting in a year with 0 stems, necessary as it would skew the model if it were not included
      bind_rows(data.frame(year = 2016, transect = 60, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "YTB")) %>%
      bind_rows(data.frame(year = 2016, transect = 61, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "YTB")) %>%
      bind_rows(data.frame(year = 2016, transect = 63, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "YTB")) %>%
      bind_rows(data.frame(year = 2017, transect = 71, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "SKY")) %>%
      bind_rows(data.frame(year = 2017, transect = 62, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "YTB")) %>%
      bind_rows(data.frame(year = 2017, transect = 65, N_seedlings = 0, N_total = 0, N_budlings = 0, herb_mean = 0, site = "YTB")) %>%
      group_by(year, transect)

    #transects 44 and 48 were abandoned after 2013 and so were not included
    data13_14 <- data_gp %>% filter(year %in% 2013:2014 & ! transect %in% c(44, 48)) %>%
      group_by(transect) %>%
      summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                herb_mean = first(herb_mean),
                site = first(site))

    #transects 70 and 72 were abandoned after 2014 and so were not included
    data14_15 <- data_gp %>% filter(year %in% 2014:2015 & ! transect %in% c(70, 72)) %>%
      group_by(transect) %>%
      summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                herb_mean = first(herb_mean),
                site = first(site))

    #transect 80 was abandoned after 2015 so was not included
    data15_16 <- data_gp %>% filter(year %in% 2015:2016 & transect != 80) %>%
      group_by(transect) %>%
      summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                herb_mean = first(herb_mean),
                site = first(site))

    #GRN was a new site and data is only available for 2017, and transect 80 in 2017 is a
    #  new transect NOT the same as transect 80 that was abandoned (same number only, which
    #  should be fixed.) Also, transects 60, 61, 62 are not present in 2017 since they died,
    #  so they are excluded.
    data16_17 <- data_gp %>% filter(year %in% 2016:2017 & site != 'GRN' & ! transect %in% c(60, 61, 63, 73, 80)) %>%
      group_by(transect) %>%
      summarize(bdlgs_per_stem = last(N_budlings)/first(N_total),
                herb_mean = first(herb_mean),
                site = first(site))

    fulldat <- bind_rows(data13_14, data14_15, data15_16, data16_17)

    merged <- fulldat %>% group_by(transect) %>%
      summarize(herb_mean = mean(herb_mean),
                bdlgs_per_stem = mean(bdlgs_per_stem),
                site = first(site))

    cat(paste0("Calculating budlings per stem fit..."))
    mdl <- lm(bdlgs_per_stem ~ herb_mean, data=merged)
    cat("done!\n")

    # Add data, as nls doesn't store the data
    mdl$merged <- merged

    budlings.per.stem.fit <- vector("list", 3)
    budlings.per.stem.fit[[1]] <- mdl
    budlings.per.stem.fit[[2]] <- merged$site
    if (obj$mdlargs$method == 'pow') {
      budlings.per.stem.fit[[3]] <- function (x, pars, perturb = rep(0,2)) {
        pars <- pars + perturb
        y <- pars[1]*x^pars[2]
      }
    } else {
      budlings.per.stem.fit[[3]] <- function (x, pars, perturb = rep(0,2)) {
        pars <- pars + perturb
        y <- pars[1] + pars[2]*x
      }
    }
    names(budlings.per.stem.fit) <- c("fit","site","predict")

    if (saveresults) {
      save(budlings.per.stem.fit, file = file.path(mwCache,"budlingsPerStemFit.RData"))
    }
  } else {
    load(file.path(mwCache,"budlingsPerStemFit.RData"))
  }

  obj$pars$budlings.per.stem.fit <- budlings.per.stem.fit

  return(obj)
}

# Matrices ------------------------------

#' Create flowering matrix.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
setFloweringMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Flowering
  # N x N^2
  obj$matrices$F = t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(x,y) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = x, herb_avg = y), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]

  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }

  return(obj)
}

#' Create survival matrix.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
setSurvivalMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Survival
  # N x N^2
  obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(x,y) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = x, herb_avg = y), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]

  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }

  return(obj)
}

#' Create growth matrix.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %>%
setGrowthMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,5)) {
  # Growth
  # N x N^2

  N <- obj$N
  Mu <- c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(x,y) {predict(obj$pars$growth.fit, newdata = data.frame(h_apical = x, herb_avg = y), type=obj$site, perturb=perturb[1:4])}))
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
  # colSums(G%*%H*dx.h_apical.next*dx.herb_avg)

  # Plot growth
  # image.plot(h_apical, h_apical.next, t(G%*%H), col=topo.colors(100))

  if (update) {
    obj <- obj %>% computeSexualKernel()
    obj <- computeMPM(obj)
  }

  return(obj)
}

#' Create pods matrix.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %<>%
setPodsMatrix.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  # Pods
  # N x N^2

  N <- obj$N
  P <- matrix(rep(0, N^3), nrow=N)
  for (i in 1:N) {
    # h_apical.next
    for (j in 1:N) {
      # herb_avg
      P[i,(1+(j-1)*N):(j*N)] <- predict(obj$pars$pods.fit, newdata = data.frame(h_apical.next = obj$vars$h_apical.next$x[i], herb_avg = obj$vars$herb_avg$x[j]), type=obj$site, perturb=perturb)
    }
  }
  obj$matrices$P <- P

  if (update) {
    obj %<>% computeSexualKernel()
    obj <- computeMPM(obj)
  }

  return(obj)
}

#' Create herbivory matrix.
#'
#' @param obj A mwIPM model object.
#' @param dist.herb Herbivory distribution.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#' @param distpars Determines in which scale (log or linear) perturbation is occurring.
#'
#' @return A mwIPM model object.
#' @export
#'
#' @importFrom magrittr %<>%
setHerbivoryMatrix.mwIPM <- function(obj, dist.herb = NA, update = TRUE, perturb = rep(0,3), distpars = FALSE) {
  # Herbivory matrix
  # N x N^2

  N <- obj$N
  if (any(is.na(dist.herb))) {
    dist.herb <- obj$pars$munched.fit[[obj$site]]$predict(obj$vars$herb_avg$b,
                                                          c(obj$pars$munched.fit[[obj$site]]$pmunch,
                                                            obj$pars$munched.fit[[obj$site]]$fit$estimate),
                                                          perturb,
                                                          distpars = distpars)
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

#' Create seedling recruitment matrix.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %<>%
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

#' Compute sexual reproduction kernel.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %<>%
computeSexualKernel.mwIPM <- function(obj, update = TRUE, perturb = rep(0,2)) {
  attach(obj$matrices, warn.conflicts = TRUE)

  Ks <- (obj$pars$seedling.emergence[[obj$site]]+perturb[1])*
    (obj$pars$seeds.per.pod+perturb[2])*R%*%(P*G*S*obj$matrices$F)%*%H*obj$vars$herb_avg$dx*obj$vars$h_apical.next$dx*obj$vars$h_apical$dx

  detach(obj$matrices)

  obj$kernels$Ks <- Ks

  if (update) {
    obj %<>% computeFullKernel()
  }

  return(obj)
}

#' Compute clonal reproduction kernel.
#'
#' @param obj A mwIPM model object.
#' @param update Update dependencies?
#' @param perturb Parameter perturbation vector for sensitivity analysis.
#'
#' @return A mwIPM model object.
#'
#' @importFrom magrittr %<>%
computeClonalKernel.mwIPM <- function(obj, update = TRUE, perturb = rep(0,4)) {
  ## Budling Recruitment (same as Kc)
  # perturb = c(budlings.per.stem, budlings)

  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)

  if (obj$mdlargs$input == 'meanonly') {
    mean.herb_avg <- sum(dist.herb*herb_avg$x)*herb_avg$dx
    mean.buds.per.stem <- budlings.per.stem.fit$predict(mean.herb_avg,
                                                        stats::coef(budlings.per.stem.fit$fit),
                                                        perturb[1:2])
  } else {

    mean.buds.per.stem <- t(budlings.per.stem.fit$predict(herb_avg$x,
                                                          stats::coef(budlings.per.stem.fit$fit),
                                                          perturb[1:2])) %*%
      t(t(H[1+(0:(obj$N-1))*obj$N,1]))*herb_avg$dx
  }

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

#' Compute full kernel.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object.
computeFullKernel.mwIPM <- function(obj) {
  obj$kernels$K <- obj$kernels$Ks + obj$kernels$Kc
  # image.plot(h_apical, h_apical, t(K), col=topo.colors(100))
  # contour(h_apical, h_apical, t(K), add = TRUE, drawlabels = TRUE)

  return(obj)
}

# Analysis ------------------------------

#' Set site and recompute kernel.
#'
#' @param obj A mwIPM model object.
#' @param site Site ("Bertha", "BLD1", "BLD2", "GET", "GRN", "PWR", "SKY", "YTB")
#' @param compute Recompute or load from cache?
#'
#' @return A mwIPM model object.
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @examples
#' ipm %<>% setSite("BLD1")
setSite.mwIPM <- function(obj, site = "Bertha", compute = FALSE) {
  if ((compute) | is.na(obj$site) | (site != obj$site)) {
    if (site %in% obj$all_sites) {
      obj$site <- site

      obj <- obj %>% setVars()

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
      stop(paste0(site, " is not a valid site!  Please choose one of ", toString(obj$all_sites), "."))
    }
  } else {
    warning(paste(site, "is already the current site."))
  }

  return(obj)
}

#' Compute the MPM from the IPM
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object.
#'
#' @export
#'
#' @examples
#' ipm %<>% computeMPM()
computeMPM.mwIPM <- function(obj) {
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)

  if (obj$mdlargs$input == 'meanonly') {
    mean.herb_avg <- sum(dist.herb*herb_avg$x)*herb_avg$dx
    mean.buds.per.stem <- budlings.per.stem.fit$predict(mean.herb_avg,
                                                        stats::coef(budlings.per.stem.fit$fit))
  } else {
    mean.buds.per.stem <- t(budlings.per.stem.fit$predict(herb_avg$x,
                                                          stats::coef(budlings.per.stem.fit$fit))) %*%
      t(t(H[1+(0:(obj$N-1))*obj$N,1]))*herb_avg$dx
  }

  Kss <- seedling.emergence[[obj$site]]*seeds.per.pod*(P*G*S*obj$matrices$F)%*%H*herb_avg$dx*h_apical.next$dx*h_apical$dx
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

#' Compute the population growth rate.
#'
#' @param obj A mwIPM model object.
#'
#' @return The population growth rate
#' @export
#'
#' @examples
#' lambda <- ipm %>% analyzeGrowthRate()
analyzeGrowthRate.mwIPM <- function(obj) {
  return(Re(eigen(obj$kernels$K)$values[1]))
}

#' Compute standard IPM population level metrics, such as the population
#' growth rate, left and right eigenvectors, and sensitivity and elasticit
#' matrices.
#'
#' @param obj A mwIPM model object.
#'
#' @return A mwIPM model object, with the results of the computation stored
#' in a list.
#' @export
#'
#' @examples
#' ipm %<>% analyzeStandard()
#' str(ipm$analysis$standard)
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

#' Sensitivity analysis on parameters.
#'
#' @param obj A mwIPM model object.
#' @param compute Recompute or load from cache?
#' @param saveresults Cache results?
#' @param params Which parameter to perturb?  One of "all", "Flowering", ...
#' @param distpars Determines in which scale (log or linear) perturbation is occurring.
#'
#' @return A mwIPM model object, with the results of the computation stored
#' in a list.
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#'
#' @examples
#' ipm %<>% analyzeParameters()
#' str(ipm$analysis$parameters)
analyzeParameters.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, params = "all", distpars = FALSE) {
  if (!file.exists(file.path(mwCache,"parameterAnalysis.RData")) | (compute)) {
    # Initialize to empty data frame
    analysis <- tbl_df(data.frame(sensitivity = NULL,
                                  pars = NULL,
                                  type = NULL,
                                  name = NULL))

    if (("all" %in% params) | ("Flowering" %in% params)) {
      cat("Analyzing flowering fit...")
      flowering_func <- function(x) {obj %>% setFloweringMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(flowering_func,rep(0,4)),
                          pars = obj$pars$flower.fit$pars$unscaled[obj$site,],
                          type = as.character("Flowering"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "h_apical:herb_avg")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Survival" %in% params)) {
      cat("Analyzing survival fit...")
      survival_func <- function(x) {obj %>% setSurvivalMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(survival_func, rep(0,4)),
                          pars = obj$pars$surv.fit$pars$unscaled[obj$site,],
                          type = as.character("Survival"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "h_apical:herb_avg")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Growth" %in% params)) {
      cat("Analyzing growth fit...")
      growth_func <- function(x) {obj %>% setGrowthMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(growth_func, rep(0,5)),
                          pars = c(obj$pars$growth.fit$pars$unscaled[obj$site,],
                                   obj$pars$growth.fit$pars$sd[obj$site]),
                          type = as.character("Growth"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "h_apical:herb_avg",
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
        tbl_df(data.frame(sensitivity = numDeriv::grad(pods_func, rep(0,4)),
                          pars = obj$pars$pods.fit$pars$unscaled[obj$site,],
                          type = as.character("Pods"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "h_apical:herb_avg")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Seedlings" %in% params)) {
      cat("Analyzing seedling distribution...")
      seedling_func <- function(x) {obj %>% setSeedlingRecruitmentMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(seedling_func, rep(0,2)),
                          pars = ParsToMoms(x = obj$pars$seedling.fit$fit$estimate,
                                            type = ifelse(!distpars,
                                                          obj$pars$seedling.fit$fit$distname,
                                                          "id")),
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
        tbl_df(data.frame(sensitivity = numDeriv::grad(sexual_func, rep(0,2)),
                          pars = c(obj$pars$seedling.emergence[obj$site],
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
        tbl_df(data.frame(sensitivity = numDeriv::grad(clonal_func, rep(0,4)),
                          pars = c(stats::coef(obj$pars$budlings.per.stem.fit$fit),
                                   obj$pars$budling.fit[[obj$site]]$fit$estimate),
                          type = c("Clonal", "Clonal", "Budlings", "Budlings"),
                          name = c("a",
                                   "b",
                                   "mean",
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Herbivory" %in% params)) {
      cat("Analyzing herbivory distribution...")
      herbivory_func <- function(x) {obj %>% setHerbivoryMatrix(perturb = x, distpars = distpars) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(herbivory_func, rep(0,3)),
                          pars = c(obj$pars$munched.fit[[obj$site]]$pmunch,
                                   ParsToMoms(x = obj$pars$munched.fit[[obj$site]]$fit$estimate,
                                              type = ifelse(!distpars,
                                                            obj$pars$munched.fit[[obj$site]]$fit$distname,
                                                            "id"))),
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
      mutate(lambda = obj %>% analyzeGrowthRate(),
             elasticity = sensitivity/(lambda/pars))
    analysis$type <- factor(analysis$type)

    if (saveresults) {
      save(analysis, file = file.path(mwCache,"parameterAnalysis.RData"))
    }
  } else {
    load(file.path(mwCache,"parameterAnalysis.RData"))
  }

  obj$analysis$parameters <- analysis

  return(obj)
}

# Renderers ------------------------------

#' Render plot of flowering vs. height.
#'
#' @param obj A mwIPM model object.
#'
#' @return A plot object.
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#' ipm %>% renderFlowerFit()
renderFloweringFit.mwIPM <- function(obj) {
  requirePackages(c("grid", "gridExtra", "gtable", "RColorBrewer"))

  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$vars, warn.conflicts = FALSE)

  ## Top (heatmap)
  M <- mesh(h_apical$x, herb_avg$x)
  plotdata <- tbl_df(data.frame(h_apical = as.vector(M$x),
                                herb_avg = as.vector(M$y)))
  z <- predict(flower.fit,
               newdata=data.frame(h_apical = as.vector(M$x),
                                  herb_avg = as.vector(M$y)),
               type="Bertha")
  plotdata <- plotdata %>% mutate(prob.flower = z)

  # Heatmap
  imgt <- ggplot(plotdata, aes(x = h_apical, y = herb_avg, z = prob.flower)) +
    geom_raster(aes(fill = prob.flower)) +
    geom_contour(colour = "white", alpha = 0.8) +
    scale_fill_gradientn("Flowering\nProbability",
                         colours=c("#00000000","#BBBBBBBB"),
                         limits=c(0, 1))

  imgt <- imgt + theme(axis.line = element_blank()) +
    scale_x_continuous(limits = c(0, h_apical$max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(herb_avg$min, herb_avg$max), expand = c(0, 0)) +
    ylab("Herbivory Score")

  # Add lines
  herb_ex <- data.frame(yintercept = c(0, mean(obj$data$herb_avg, na.rm=T), 6))

  imgt <- imgt + geom_hline(aes(yintercept = yintercept),
                          data = herb_ex,
                          linetype = c(1, 2, 5),
                          col = RColorBrewer::brewer.pal(5, "YlOrRd")[2:4],
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
                                                             herb_avg = herb_ex[1,1]),
                                          type="Bertha"),
                    Herbivory = "min")
  avg <- data.frame(h_apical = h_apical$x,
                    prob.flower = predict(flower.fit,
                                          newdata=data.frame(h_apical = h_apical$x,
                                                             herb_avg = herb_ex[2,1]),
                                          type="Bertha"),
                    Herbivory = "avg")
  max <- data.frame(h_apical = h_apical$x,
                    prob.flower = predict(flower.fit,
                                          newdata=data.frame(h_apical = h_apical$x,
                                                             herb_avg = herb_ex[3,1]),
                                          type="Bertha"),
                    Herbivory = "max")

  mycurves <- rbind(min, avg, max)
  mycurves$Herbivory <- as.factor(mycurves$Herbivory)

  imgb <- ggplot(mycurves, aes(x = h_apical, y = prob.flower, col = Herbivory, linetype = Herbivory)) +
    geom_line(size = 1.0) +
    scale_x_continuous(limits = c(0, h_apical$max), expand = c(0, 0)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlOrRd")[2:4]) +
    xlab("Apical Height (cm)") +
    ylab("Flowering Probability")

  imgb <- imgb + theme(legend.position = "none")

  gt <- ggplotGrob(imgt)
  gb <- ggplotGrob(imgb)
  gb <- gtable::gtable_add_cols(gb, unit(1, "mm"))
  g <- rbind(gt, gb, size="first")
  g$widths <- grid::unit.pmax(gt$widths, gb$widths)

  img <- gridExtra::grid.arrange(g)

  detach(obj$pars)
  detach(obj$vars$h_apical)

  return(img)
}

#' Render plot of budling distributions.
#'
#' @param obj A mwIPM model object.
#'
#' @return A plot object.
#' @export
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' ipm %>% renderBudlingDistFit()
renderBudlingDistFit.mwIPM <- function(obj) {
  requirePackages("scales")

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
    scale_fill_manual(values = c(scales::hue_pal()(5), NA)) +
    scale_color_manual(values = c(scales::hue_pal()(5), "black")) +
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

#' Render plot of herbivory distributions.
#'
#' @param obj A mwIPM model object.
#'
#' @return A plot object.
#' @export
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' ipm %>% renderHerbivoryDistFit()
renderHerbivoryDistFit.mwIPM <- function(obj) {
  requirePackages("scales")

  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$vars$herb_avg, warn.conflicts = FALSE)

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
    xlab("Herbivory Score") +
    ylab("Probability Density") +
    scale_x_continuous(limits = c(0, 6)) +
    scale_fill_manual(values = c(scales::hue_pal()(5), NA)) +
    scale_color_manual(values = c(scales::hue_pal()(5), "black")) +
    labs(colour="Sites", fill="Sites") +
    theme(legend.background = element_rect(fill="lightgrey",
                                           size=0.1,
                                           linetype="solid"),
          legend.key.size =  unit(0.18, "in"),
          legend.position = c(0.885, 0.72))

  detach(obj$pars)
  detach(obj$vars$herb_avg)

  return(pa)
}

# Helpers  ------------------------------

#' Parameters to moments.
#'
#' @param x Vector of parameters.
#' @param type Distribution type (lnorm, gamma, or identity)
#'
#' @return Vector of moments.
ParsToMoms <- function(x, type = "ident") {
  if (type == "lnorm") {
    y <- c(y1 = exp(x[1] + 0.5*x[2]*x[2]),
           y2 = sqrt(exp(2*x[1] + x[2]*x[2])*(exp(x[2]*x[2])-1)))
    names(y) <- c("mean","sd")
  } else if (type == "gamma") {
    y <- c(y1 = x[1]/x[2],
           y2 = sqrt(x[1])/x[2])
    names(y) <- c("mean","sd")
  } else {
    y <- x
    names(y) <- names(x)
  }
  return(y)
}

#' Jacobian of parameters to moments.
#'
#' @param x Vector of parameters.
#' @param type Distribution type (lnorm, gamma, or identity)
#'
#' @return Jacobian matrix
jacParsToMoms <- function(x, type = "ident") {
  if (type == "lnorm") {
    CV <- x[2]/x[1]
    a <- 1 + CV*CV
    b <- sqrt(log(a))
    J <- (1/(x[1]*a))*diag(c(1,CV/b))%*%matrix(c(a+CV*CV, -1*CV, -1*CV, 1), byrow=T, nrow=2)
    rownames(J) <- c("mu", "sd")
    colnames(J) <- c("m", "s")
  } else {
    J <- diag(2)
    rownames(J) <- c("mu", "sd")
    colnames(J) <- c("mu", "sd")
  }
  return(J)
}

#' Moments to parameters.
#'
#' @param y Vector of moments.
#' @param type Distribution type (lnorm, gamma, or identity)
#'
#' @return Vector of moments.
MomsToPars <- function(y, type = "ident") {
  if (type == "lnorm") {
    x <- c(x1 = log(y[1]/sqrt(1 + (y[2]/y[1])^2)),
           x2 = sqrt(log(1 + (y[2]/y[1])^2)))
    names(x) <- c("meanlog","sdlog")
  } else if (type == "gamma") {
    x <- c(x1 = (y[1]/y[2])^2,
           x2 = y[1]/(y[2]*y[2]))
    names(x) <- c("shape","rate")
  } else {
    x <- y
    names(x) <- names(y)
  }
  return(x)
}

#' Jacobian of moments to parameters.
#'
#' @param y Vector of moments.
#' @param type Distribution type (lnorm, gamma, or identity)
#'
#' @return Jacobian matrix
jacMomsToPars <- function(y, type = "ident") {
  if (type == "lnorm") {
    J <- exp(y[1] + y[2]*y[2]/2)*diag(c(1,sqrt(exp(y[2]*y[2])-1)))%*%matrix(c(1,1,1,1+1/(1-exp(-1*y[2]*y[2]))), byrow=T, nrow=2)%*%diag(c(1,y[2]))
    rownames(J) <- c("m", "s")
    colnames(J) <- c("mu", "sd")
  } else {
    J <- diag(2)
    rownames(J) <- c("m", "s")
    colnames(J) <- c("m", "s")
  }
  return(J)
}

#' Perturbation of the parameter-to-moments transformation.
#'
#' @param pars Vector of parameters or moments
#' @param perturb Perturbation vector
#' @param type Distribution type (lnorm, gamma, or identity)
#'
#' @return Gradient vector.
perturbTrans <- function(pars, perturb = rep(0,2), type = "ident") {
  if (type == "lnorm") {
    MSD <- ParsToMoms(x = pars, type)
    tpars <- MomsToPars(y = MSD+perturb, type)
  } else if (type == "gamma") {
    MSD <- ParsToMoms(x = pars, type)
    tpars <- MomsToPars(y = MSD+perturb, type)
  } else {
    tpars <- pars + perturb
  }
  return(tpars)
}

#' Helper function to check for required packages.
#'
#' @param req Vector of required packages.
#'
#' @export
requirePackages <- function(req) {
  lapply(req, function(pkg) {
    if (!requireNamespace(pkg)) {
      stop(paste("Package", pkg, "is required.  Please install it to continue."), call. = FALSE)
    }
  })

  invisible()
}
