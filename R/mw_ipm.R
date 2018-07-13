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
setPodsFit <- function(obj, compute, saveresults, update) UseMethod("setPodsFit")

## Distributions
setSeedlingDistFit <- function(obj, compute, saveresults, update) UseMethod("setSeedlingDistFit")
setBudlingDistFit <- function(obj, compute, saveresults, update) UseMethod("setBudlingDistFit")
setHerbivoryDistFit <- function(obj, compute, saveresults, update) UseMethod("setHerbivoryDistFit")
setCardenolideDistFit <- function(obj, compute, saveresults, update) UseMethod("setCardenolideDistFit")
setLMADistFit <- function(obj, compute, saveresults, update) UseMethod("setLMADistFit")

## Regressions
setBudlingsPerStemFit <- function(obj, compute, saveresults, update) UseMethod("setBudlingsPerStemFit")

# Set matrices from fits
setFloweringMatrix <- function(obj, update, perturb) UseMethod("setFloweringMatrix")
setSurvivalMatrix <- function(obj, update, perturb) UseMethod("setSurvivalMatrix")
setPodsMatrix <- function(obj, update, perturb) UseMethod("setPodsMatrix")
setVarDistMatrices <- function(obj, update, perturb) UseMethod("setVarDistMatrices")
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
#mwCache <- file.path(rappdirs::app_dir("Milkweed", "LaMar")$data(), "calculated") # Commented out when new Cache created explicitly, same with below line
#options(milkweed.cache = mwCache) # Store this for outside use in vignettes

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
    x$data = list(stem = stemdata %>% filter(site %in% x$all_sites),
                  seed = seeddata,
                  dem = dem.data,
                  trait = trait.data,
                  full = full.data,
                  map = map.data,
                  bud = bud.data)
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
    saveresults = (length(list.files('data/mwCache', pattern=".RData")) == 0)
  } else if (!is.logical(x$saveresults)) {
    stop(paste("Variable 'saveresults' must be TRUE or FALSE."))
  } else {
    saveresults <- x$saveresults
  }
  if (all(names(x) != "mdlargs")) {
    x$mdlargs <- list(method = 'pow',
                      input = 'full')
  } else {
    if (all(names(x$mdlargs) != "method")) {
      x$mdlargs$method = "pow"
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
                                  herb_avg = list(min = NA,
                                                  max = NA,
                                                  b = NA,
                                                  x = NA,
                                                  dx = NA),
                                  LMA = list(min = NA,
                                             max = NA,
                                             b = NA,
                                             x = NA,
                                             dx = NA),
                                  Cardenolides = list(min = NA,
                                                      max = NA,
                                                      b = NA,
                                                      x = NA,
                                                      dx = NA)),
                      pars = list(flower.fit = NA,     # vital rate
                                  surv.fit = NA,       # vital rate
                                  pods.fit = NA,       # vital rate
                                  seedling.fit = NA,   # height dist
                                  budling.fit = NA,    # height dist
                                  munched.fit = NA,    # herbivory dist
                                  lma.fit = NA,        # trait dist
                                  card.fit = NA,       # trait dist
                                  budlings.per.stem.fit = NA, # vital rate
                                  seedling.emergence = NA,    # constant (actually recruitment of establishment in sexual pathway)
                                  seeds.per.pod = NA),         # constant
                      matrices = list(F = NA,
                                      S = NA,                 # F, S, and P are the vital rate matrices, N x N^4
                                      P = NA,
                                      H = NA,
                                      C = NA,                 # H, C, and L are the variable distribution matrices, N^4 x N
                                      L = NA,
                                      R = NA),                # Seedling Recruitment Matrix, N x N
                      kernels = list(Ks = NA,
                                     Kc = NA,
                                     K = NA),
                      MPM = NA,
                      analysis = list(standard = NA,
                                      parameters = NA)))

  y <- structure(x, class = "mwIPM") %>%
    setPars(compute = compute, saveresults = saveresults, update = FALSE) %>%
    setSite(compute = compute, site = site) #%>% #commented out by Soren for testing
    #computeMPM()

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
  obj$data <- list(stem = sample_frac(obj$data_orig$stem, replace=T),
                   seed = sample_frac(obj$data_orig$seed, replace=T),
                   dem = sample_frac(obj$data_orig$dem, replace=T),
                   trait = obj$data_orig$trait,
                   full = sample_frac(obj$data_orig$full, replace=T),
                   map = obj$data_orig$map,
                   bud = sample_frac(obj$data_orig$bud, replace=T))
  obj %<>% setPars(compute = TRUE, update = FALSE) %>% # why FALSE?
           setSite(compute = TRUE) #%>%
           #computeMPM() ## commented out by Soren

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
    data <- obj$data$dem %>% filter(site == obj$site)
  }
  data <- obj$data$dem %>% summarize(max_h_apical = max(h_apical, na.rm=T))

  # h_apical
  h_apical <- list(min = 0,
                   max = 1.1*data$max_h_apical)
  h_apical$b = h_apical$min + c(0:N)*(h_apical$max - h_apical$min)/N # boundary points
  h_apical$x = 0.5*(h_apical$b[1:N]+h_apical$b[2:(N+1)]) # midpoints of classes
  h_apical$dx = h_apical$b[2]-h_apical$b[1] # class size

  # herb_avg
  herb_avg <- list(min = 0.0,
                   max = 6.0)
  herb_avg$b = herb_avg$min + c(0:N)*(herb_avg$max - herb_avg$min)/N # boundary points
  herb_avg$x = 0.5*(herb_avg$b[1:N] + herb_avg$b[2:(N+1)]) # midpoints of classes
  herb_avg$dx = herb_avg$b[2] - herb_avg$b[1] # class size

  data <- obj$data$full %>% summarize(min_Card = min(Cardenolides, na.rm=T),
                                      max_Card = max(Cardenolides, na.rm=T),
                                      min_LMA = min(LMA, na.rm=T),
                                      max_LMA = max(LMA, na.rm=T))

  # Cardenolides
  Cardenolides <- list(min = 0.9*data$min_Card,
                       max = 1.1*data$max_Card)
  Cardenolides$b = Cardenolides$min + c(0:N)*(Cardenolides$max - Cardenolides$min)/N # boundary points
  Cardenolides$x = 0.5*(Cardenolides$b[1:N] + Cardenolides$b[2:(N+1)]) # midpoints of classes
  Cardenolides$dx = Cardenolides$b[2] - Cardenolides$b[1] # class size

  # LMA
  LMA <- list(min = 0.9*data$min_LMA,
              max = 1.1*data$max_LMA)
  LMA$b = LMA$min + c(0:N)*(LMA$max - LMA$min)/N # boundary points
  LMA$x = 0.5*(LMA$b[1:N] + LMA$b[2:(N+1)]) # midpoints of classes
  LMA$dx = LMA$b[2] - LMA$b[1] # class size

  obj$vars = list(h_apical = h_apical,
                  herb_avg = herb_avg,
                  LMA = LMA,
                  Cardenolides = Cardenolides)

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
                                     update = FALSE) %>%
                 setCardenolideDistFit(compute = compute,
                                       saveresults = saveresults,
                                       update = FALSE) %>%
                 setLMADistFit(compute = compute,
                               saveresults = saveresults,
                               update = FALSE)

  if (update) {
    obj <- obj %>% setFloweringMatrix() %>%
                   setSurvivalMatrix() %>%
                   setPodsMatrix() %>%
                   setVarDistMatrices(update=FALSE) %>% # WHY FALSE?
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
  if (!file.exists(file.path('data/mwCache',"seedsPerPodConst.RData")) | (compute)) {
    seeds_per_pod_data <- obj$data$seed

    cat("Computing seeds per pod constant...")
    seeds.per.pod = mean(seeds_per_pod_data$total_seed)
    cat("done!\n")

    if (saveresults) {
      save(seeds.per.pod, file = file.path('data/mwCache',"seedsPerPodConst.RData"))
    }
  } else {
    load(file.path('data/mwCache',"seedsPerPodConst.RData"))
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
  if (!file.exists(file.path('data/mwCache',"seedlingEmergenceConst.RData")) | (compute)) {
    seedling.emergence <- NA

    seeds.per.pod <- obj$pars$seeds.per.pod

    # Bertha
    cat("Computing seedling emergence constants...")
    data_gp <- obj$data$stem %>% group_by(year) %>% summarize(N_seedlings = sum(seedling, na.rm=T),
                                                              N_seeds = seeds.per.pod*sum(N_pods, na.rm=T))

    seedling.emergence <- mean(data_gp$N_seedlings[1:nrow(data_gp)-1]/data_gp$N_seeds[2:nrow(data_gp)])
    cat("done!\n")

    if (saveresults) {
      save(seedling.emergence, file = file.path('data/mwCache',"seedlingEmergenceConst.RData"))
    }
  } else {
    load(file.path('data/mwCache',"seedlingEmergenceConst.RData"))
  }

  obj$pars$seedling.emergence <- seedling.emergence

  return(obj)
}

# Fits ------------------------------

## Vital rates

#' Sets the flowering model.
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
  if (!file.exists(file.path('data/mwCache',"flowerFit.RData")) | (compute)) {
    metadata_usc <- obj$data$full %>% filter(!is.na(h_apical),
                                        !is.na(herb_avg),
                                        !is.na(fec.flower),
                                        !is.na(Cardenolides)) %>%
      select(-c(C, ADF, ADL, Cellulose, NDWI, PRI, Chl_g_m2))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, herb_avg, Cardenolides),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing flowering fit...")
    flower.mdl <- glm(fec.flower ~ h_apical + herb_avg + Cardenolides, data = metadata_sc, family = binomial())
    cat("done!\n")

    flower.fit <- mwMod(list(mdl = flower.mdl,
                             vars = c("h_apical", "herb_avg", "Cardenolides"),
                             scaled = list(h_apical = scale(metadata_usc$h_apical),
                                           herb_avg = scale(metadata_usc$herb_avg),
                                           Cardenolides = scale(metadata_usc$Cardenolides))))

    # Check parameters
    cat("Checking parameters:\n")
    checkPars(flower.fit)

    if (saveresults) {
      save(flower.fit, file = file.path('data/mwCache',"flowerFit.RData"))
    }
  } else {
    load(file.path('data/mwCache',"flowerFit.RData"))
  }

  obj$pars$flower.fit <- flower.fit

  return(obj)
}

#' Sets the survival model.
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
  if (!file.exists(file.path('data/mwCache',"survivalFit.RData")) | (compute)) {
    metadata_usc <- obj$data$full %>% filter(!is.na(h_apical),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        !is.na(surv),
                                        !is.na(LMA),
                                        !is.na(Cardenolides)) %>%
      select(-c(C, ADF, ADL, Cellulose, NDWI, PRI, Chl_g_m2))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, herb_avg, Cardenolides, LMA),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing survival fit...")
    surv.mdl <- glm(surv ~ h_apical + herb_avg + Cardenolides + LMA, data = metadata_sc, family = binomial())
    cat("done!\n")

    surv.fit <- mwMod(list(mdl = surv.mdl,
                           vars = c("h_apical", "herb_avg", 'Cardenolides', 'LMA'),
                           scaled = list(h_apical = scale(metadata_usc$h_apical),
                                         herb_avg = scale(metadata_usc$herb_avg),
                                         Cardenolides = scale(metadata_usc$Cardenolides),
                                         LMA = scale(metadata_usc$LMA)),
                           perturb.mdl = surv.mdl))

    #test code to see if manual perturbation of beta_c in survival vr function changes lambda how we think it will
    #obj$pars$surv.fit$perturb.mdl$coefficients[4] = obj$pars$surv.fit$mdl$coefficients[4] + 16

    # Check parameters
    cat("Checking parameters:\n")
    checkPars(surv.fit)

    if (saveresults) {
      save(surv.fit, file = file.path('data/mwCache',"survivalFit.RData"))
    }
  } else {
    load(file.path('data/mwCache',"survivalFit.RData"))
  }

  obj$pars$surv.fit <- surv.fit

  return(obj)
}

#' Sets the pods model.
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
  if (!file.exists(file.path('data/mwCache',"podsFit.RData")) | (compute)) {
    # We only want stems that flowered, survived and have data for pods
    metadata_usc <- obj$data$full %>% filter(!is.na(h_apical),
                                        !is.na(h_apical.next),
                                        !is.na(herb_avg),
                                        fec.flower == 1,
                                        surv == 1,
                                        !is.na(N_pods),
                                        !is.na(Cardenolides),
                                        !is.na(LMA))

    metadata_sc <- metadata_usc %>% mutate_at(.vars = vars(h_apical, herb_avg, Cardenolides, LMA),
                                              .funs = funs(as.numeric(scale(.))))

    cat("Computing pods fit...")
    pods.mdl <- glm(N_pods ~ h_apical + herb_avg + LMA + Cardenolides, data = metadata_sc, family = poisson())
    cat("done!\n")

    pods.fit <- mwMod(list(mdl = pods.mdl,
                           vars = c("h_apical", "herb_avg", "Cardenolides", "LMA"),
                           scaled = list(h_apical = scale(metadata_usc$h_apical),
                                         herb_avg = scale(metadata_usc$herb_avg),
                                         Cardenolides = scale(metadata_usc$Cardenolides),
                                         LMA = scale(metadata_usc$LMA)),
                           perturb.mdl = pods.mdl))
    # Check parameters
    cat("Checking parameters:\n")
    checkPars(pods.fit)

    if (saveresults) {
      save(pods.fit, file = file.path('data/mwCache',"podsFit.RData"))
    }
  } else {
    load(file.path('data/mwCache',"podsFit.RData"))
  }

  obj$pars$pods.fit <- pods.fit

  return(obj)
}

## Distributions

#' Sets the seedling distribution.
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
  if (!file.exists(file.path('data/mwCache',"seedlingDistFit.RData")) | (compute)) {
    h_apical <- (obj$data$dem %>% filter(seedling == 1, !is.na(h_apical)))$h_apical

    cat("Computing seedling distribution fit...")
    f0 <- fitdistrplus::fitdist(h_apical, "norm")
    cat("done!\n")

    seedling.fit <- vector("list", 2)
    seedling.fit[[1]] <- f0

    seedling.fit[[2]] <- function(x, pars, perturb = rep(0,2)) {
           pars = pars + perturb
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
      save(seedling.fit, file = file.path('data/mwCache',"seedlingDistFit.RData"))
    }

  } else {
    load(file.path('data/mwCache',"seedlingDistFit.RData"))
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
  if (!file.exists(file.path('data/mwCache',"budlingDistFit.RData")) | (compute)) {
    h_apical <- (obj$data$dem %>% filter(seedling == 0,
                                     !is.na(h_apical)))$h_apical

    cat("Computing budling distribution fit...")
    f0 <- fitdistrplus::fitdist(h_apical, 'norm')
    cat("done!\n")

    budling.fit <- vector("list", 2)
    budling.fit[[1]] <- f0
    budling.fit[[2]] <- function(x, pars, perturb = rep(0,2)) {
      pars = pars + perturb
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
    names(budling.fit) <- c("fit", "predict")

    if (saveresults) {
      save(budling.fit, file = file.path('data/mwCache',"budlingDistFit.RData"))
    }

  } else {
    load(file.path('data/mwCache',"budlingDistFit.RData"))
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
  if (!file.exists(file.path('data/mwCache',"herbivoryDistFit.RData")) | (compute)) {
    metadata = obj$data$dem %>% filter(!is.na(h_apical),
                                       !is.na(munched))
    herb_avg <- (metadata %>% filter(munched == 1))$herb_avg
    pmunch <- sum(metadata$munched == 1)/nrow(metadata)

    cat("Computing herbivory distribution fit for...")
    f1 <- fitdistrplus::fitdist(herb_avg, "lnorm")
    cat('done!\n')

    munched.fit <- vector('list', 3)
    munched.fit[[1]] <- f1
    munched.fit[[2]] <- pmunch
    munched.fit[[3]] <- function(x, pars, perturb = rep(0,3), justmunch=FALSE) {
      pars[1] = pars[1] + perturb[1]
      # pars <- perturbTrans(pars, perturb)  #two lines below replace perturbTrans() call
      MSD <- ParsToMoms(x = pars[2:3], type = "lnorm")
      pars[2:3] <- MomsToPars(y = MSD+perturb[2:3], type = "lnorm")
      N <- length(x)
      dx <- x[2]-x[1]
      y <- rep(0, N-1)
      for (j in 1:(N-1)) {
        y[j] = plnorm(x[j+1], pars[2], pars[3]) - plnorm(x[j], pars[2], pars[3])
      }
      y[1] <- y[1] + plnorm(x[1], pars[2], pars[3])
      y[N-1] <- y[N-1] + plnorm(x[N], pars[2], pars[3], lower.tail=FALSE)
      y <- pars[1]*y
      if (!justmunch) {
        y[1] <- y[1] + (1-pars[1])
      }
      y <- y/dx
    }
    names(munched.fit) <- c("fit", "pmunch", "predict")

    if (saveresults) {
      save(munched.fit, file = file.path('data/mwCache',"herbivoryDistFit.RData"))
    }

  } else {
    load(file.path('data/mwCache',"herbivoryDistFit.RData"))
  }

  obj$pars$munched.fit <- munched.fit

  return(obj)
}

#' Sets the cardenolide distribution.
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
setCardenolideDistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path('data/mwCache',"cardenolideDistFit.RData")) | (compute)) {
    metadata = obj$data$full %>% filter(!is.na(Cardenolides))

    Cardenolides <- metadata$Cardenolides

    cat("Computing cardenolide distribution fit for...")
    f0 <- fitdistrplus::fitdist(Cardenolides, "norm")
    cat('done!\n')

    card.fit <- vector('list', 2)
    card.fit[[1]] <- f0
    card.fit[[2]]  <- function(x, pars, perturb = rep(0,2)) {
      pars <- pars + perturb
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
    names(card.fit) <- c("fit", "predict")

    if (saveresults) {
      save(card.fit, file = file.path('data/mwCache',"cardenolideDistFit.RData"))
    }

  } else {
    load(file.path('data/mwCache',"cardenolideDistFit.RData"))
  }

  obj$pars$card.fit <- card.fit

  return(obj)
}

#' Sets the LMA distribution.
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
setLMADistFit.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE) {
  if (!file.exists(file.path('data/mwCache',"lmaDistFit.RData")) | (compute)) {
    metadata = obj$data$full %>% filter(!is.na(LMA))

    LMA <- metadata$LMA

    cat("Computing LMA distribution fit for...")
    f0 <- fitdistrplus::fitdist(LMA, "norm")
    cat('done!\n')

    lma.fit <- vector('list', 2)
    lma.fit[[1]] <- f0
    lma.fit[[2]]  <- function(x, pars, perturb = rep(0,2)) {
      pars <- pars + perturb
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
    names(lma.fit) <- c("fit", "predict")

    if (saveresults) {
      save(lma.fit, file = file.path('data/mwCache',"lmaDistFit.RData"))
    }

  } else {
    load(file.path('data/mwCache',"lmaDistFit.RData"))
  }

  obj$pars$lma.fit <- lma.fit

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
  if (!file.exists(file.path('data/mwCache',"budlingsPerStemFit.RData")) | (compute)) {

    ## Results from Transect-Blocking.Rmd
    bud.data = obj$data$bud %>% filter(!is.na(buds_per_stem),
                                       !is.na(herb_mean),
                                       !is.na(Card_mean)) %>% #all the traits have NAs in identical rows)
                                mutate(Cardenolides = Card_mean) #makes new column to keep names the same for predictions in clonal kernal

    cat(paste0("Calculating budlings per stem fit..."))
    mdl = lm(log(buds_per_stem) ~ Cardenolides, data=bud.data)
    cat("done!\n")

    budlings.per.stem.fit <- vector("list", 2)
    budlings.per.stem.fit[[1]] = mdl
    budlings.per.stem.fit[[2]] <- function (x, pars, perturb = rep(0,2)) {
      pars <- pars + perturb
      y <- exp(pars[1] + pars[2]*x) #reversing log from lm
    }

    names(budlings.per.stem.fit) <- c("fit","predict")

    # # Using nls instead of lm even if linear for cleaner code - fits are the same
    # cat(paste0("Calculating budlings per stem fit (", obj$mdlargs$method, ")..."))
    # if (obj$mdlargs$method == 'pow') {
    #   mdl <- nls(bdlgs_per_stem ~ a*herb_mean^b, data=merged, start=list(a = 1, b = 1))
    # } else {
    #   mdl <- nls(bdlgs_per_stem ~ a + b*herb_mean, data=merged, start=list(a = 1, b = 1))
    # }
    # cat("done!\n")
    #
    # # Add data, as nls doesn't store the data
    # mdl$merged <- merged
    #
    # budlings.per.stem.fit <- vector("list", 3)
    # budlings.per.stem.fit[[1]] <- mdl
    # budlings.per.stem.fit[[2]] <- merged$site
    # budlings.per.stem.fit[[3]] <- function (x, pars, perturb = rep(0,2)) {
    #     pars <- pars + perturb
    #     y <- pars[1]*exp(pars[2])*x
    #   }
    #
    # names(budlings.per.stem.fit) <- c("fit","site","predict")

    if (saveresults) {
      save(budlings.per.stem.fit, file = file.path('data/mwCache',"budlingsPerStemFit.RData"))
    }
  } else {
    load(file.path('data/mwCache',"budlingsPerStemFit.RData"))
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
setFloweringMatrix.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE, perturb = rep(0,4)) {
  # Flowering
  # N x N^4

  if (!file.exists(file.path('data/mwCache',"flowerMat.RData")) | (compute)) {

    cat(paste0("Computing flowering matrix..."))
    cards = obj$vars$Cardenolides$x
    #start off the matrix for the first card class
    obj$matrices$F = t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(h,w) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = h, herb_avg = w, Cardenolides = cards[1]), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]
    #then recompute the y values for each card class and append to the matrix horizontally
    for (i in 2:obj$N) {
      obj$matrices$F = cbind(obj$matrices$F, t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(h,w) {predict(obj$pars$flower.fit, newdata = data.frame(h_apical = h, herb_avg = w, Cardenolides = cards[i]), type=obj$site, perturb=perturb)})))[rep(1,obj$N),])
    }
    #need to repeat the whole thing from above N times horizontally for LMA classes
    obj$matrices$F = do.call("cbind", rep(list(obj$matrices$F), obj$N))
    cat(paste0("done!\n"))

    if (update) {
      obj <- obj %>% computeSexualKernel()
      #obj <- computeMPM(obj)
    }

    if (saveresults) {
      save(obj$matrices$F, file = file.path('data/mwCache',"flowerMat.RData"))
    }

  } else {
    load(file.path('data/mwCache',"flowerMat.RData"))
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
setSurvivalMatrix.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE, perturb = rep(0,5)) {
  # Survival
  # N x N^4

  N = obj$N
  n = N^4

  if (!file.exists(file.path('data/mwCache',"survivalMat.RData")) | (compute)) {

    cat(paste0("Computing survival matrix..."))
    # lma = obj$vars$LMA$x
    # cards = obj$vars$Cardenolides$x
    # #start off the matrix for the first card and herb class
    # obj$matrices$S = t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(h,w) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = h, herb_avg = w, Cardenolides = cards[1], LMA = lma[1]), type=obj$site, perturb=perturb)})))[rep(1,obj$N),]
    # for (i in 2:obj$N) {
    #   obj$matrices$S = cbind(obj$matrices$S, t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(h,w) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = h, herb_avg = w, Cardenolides = cards[i], LMA = lma[1]), type=obj$site, perturb=perturb)})))[rep(1,obj$N),])
    # }
    # #then recompute y values for each card and lma class and append them to the matrix horizontally
    # for (j in 2:obj$N) {
    #   cat(paste0("j = ", j, "\n")) #testing progress
    #   for (i in 2:obj$N) {
    #     obj$matrices$S = cbind(obj$matrices$S, t(c(outer(obj$vars$h_apical$x, obj$vars$herb_avg$x, function(h,w) {predict(obj$pars$surv.fit, newdata = data.frame(h_apical = h, herb_avg = w, Cardenolides = cards[i], LMA = lma[j]), type=obj$site, perturb=perturb)})))[rep(1,obj$N),])
    #   }
    # }

    # combinations, going from h being first inner-most (left) layer to l being last outer-most (right) layer
    mtx <- data.matrix(
      expand.grid(
        h = obj$vars$h_apical$x,
        w = obj$vars$herb_avg$x,
        c = obj$vars$Cardenolides$x,
        l = obj$vars$LMA$x
      )
    )

    # using scaled version for the faster predict method
    hs <- (obj$vars$h_apical$x - attr(obj$pars$surv.fit$scaled$h_apical,"scaled:center"))/attr(obj$pars$surv.fit$scaled$h_apical,"scaled:scale")
    ws <- (obj$vars$herb_avg$x - attr(obj$pars$surv.fit$scaled$herb_avg,"scaled:center"))/attr(obj$pars$surv.fit$scaled$herb_avg,"scaled:scale")
    cs <- (obj$vars$Cardenolides$x - attr(obj$pars$surv.fit$scaled$Cardenolides,"scaled:center"))/attr(obj$pars$surv.fit$scaled$Cardenolides,"scaled:scale")
    ls <- (obj$vars$LMA$x - attr(obj$pars$surv.fit$scaled$LMA,"scaled:center"))/attr(obj$pars$surv.fit$scaled$LMA,"scaled:scale")
    mtx <- data.matrix(expand.grid(h_apical = hs, herb_avg = ws, Cardenolides = cs, LMA = ls))

    #copy model, coefficients + perturb --> into predict
    obj$pars$surv.fit$perturb.mdl$coefficients = obj$pars$surv.fit$perturb.mdl$coefficients + perturb
    K <- predict(obj$pars$surv.fit$perturb.mdl, newdata = data.frame(mtx), type="response")

    S = t(K)[rep(1,N),]
    obj$matrices$S <- S
    cat(paste0("done!\n"))

    if (update) {
      obj <- obj %>% computeSexualKernel()
      #obj <- computeMPM(obj)
    }

    if (saveresults) {
      save(obj$matrices$S, file = file.path('data/mwCache',"survivalMat.RData"))
    }

  } else {
    load(file.path('data/mwCache',"survivalMat.RData"))
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
setPodsMatrix.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE, perturb = rep(0,5)) {
  # Pods
  # N x N^4

  N = obj$N
  n = N^4

  if (!file.exists(file.path('data/mwCache',"podsMat.RData")) | (compute)) {

    cat(paste0("Computing pods matrix..."))
    # combinations, going from h being first inner-most (left) layer to l being last outer-most (right) layer
    mtx <- data.matrix(
      expand.grid(
        h = obj$vars$h_apical$x,
        w = obj$vars$herb_avg$x,
        c = obj$vars$Cardenolides$x,
        l = obj$vars$LMA$x
      )
    )

    # Scaled version
    hs <- (obj$vars$h_apical$x - attr(obj$pars$pods.fit$scaled$h_apical,"scaled:center"))/attr(obj$pars$pods.fit$scaled$h_apical,"scaled:scale")
    ws <- (obj$vars$herb_avg$x - attr(obj$pars$pods.fit$scaled$herb_avg,"scaled:center"))/attr(obj$pars$pods.fit$scaled$herb_avg,"scaled:scale")
    cs <- (obj$vars$Cardenolides$x - attr(obj$pars$pods.fit$scaled$Cardenolides,"scaled:center"))/attr(obj$pars$pods.fit$scaled$Cardenolides,"scaled:scale")
    ls <- (obj$vars$LMA$x - attr(obj$pars$pods.fit$scaled$LMA,"scaled:center"))/attr(obj$pars$pods.fit$scaled$LMA,"scaled:scale")
    mtx <- data.matrix(expand.grid(h_apical = hs, herb_avg = ws, Cardenolides = cs, LMA = ls))

    #copy model, coefficients + perturb --> into predict
    obj$pars$pods.fit$perturb.mdl$coefficients = obj$pars$pods.fit$perturb.mdl$coefficients + perturb
    K <- predict(obj$pars$pods.fit$perturb.mdl, newdata = data.frame(mtx), type="response", re.form=NA)

    P = t(K)[rep(1,N),]

    # K <- array(0, dim=c(N,n))  # initialize with zeros

    # imtx <- data.matrix( expand.grid(h=1:N, w=1:N, c=1:N, l=1:N) )  # Matrix of indices

    #K[imtx] <- apply(mtx, 1, function(x) {
    #  predict(obj$pars$pods.fit, newdata = data.frame(h_apical = x["h"], herb_avg = x["w"], Cardenolides = x["c"], LMA = x["l"]))
    # predict(obj$pars$pods.fit, newdata = data.frame(h_apical = x["h"], herb_avg = x["w"], card = x["c"], lma = x["l"]))
    # cat(paste0('h_apical = ', x["h"], ' herb_avg = ', x["w"], ' card = ', x["c"], ' lma = ', x["l"], '\n', ' y = ', y)) #test output
    #  })
    # c(K) will give you ONE row of the matrix that you want.  Replicate this row N times
    # P = c(K)[rep(1,N),]

    obj$matrices$P <- P
    cat(paste0("done!\n"))

    if (update) {
      obj %<>% computeSexualKernel()
      #obj <- computeMPM(obj)
    }

    if (saveresults) {
      save(obj$matrices$P, file = file.path('data/mwCache',"podsMat.RData"))
    }

  } else {
    load(file.path('data/mwCache',"podsMat.RData"))
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
setVarDistMatrices.mwIPM <- function(obj, compute = FALSE, saveresults = FALSE, update = TRUE, perturb = rep(0,7)) {
  # Herbivory matrix
  # N^4 x N

  if (!file.exists(file.path('data/mwCache',"varDistMats.RData")) | (compute)) {

    library(pracma)
    attach(obj$vars, warn.conflicts = FALSE)
    attach(obj$pars, warn.conflicts = FALSE)

    N = obj$N
    n = N^4
    iA <- cbind(1:n, 1+pracma::mod(0:(n-1), N)) #matrix of indices

    #combinations, going from h being first inner-most (left) layer to l being last outer-most (right) layer
    cat('expanding grid and initializing...')
    mtx <- data.matrix(
      expand.grid(
        h = h_apical$x, #the value of this doesn't actually matter, we just need the dimensions to match properly (need combinations of height)... see ?expand.grid
        w = munched.fit$predict(herb_avg$b,
                                c(munched.fit$pmunch,
                                  munched.fit$fit$estimate),
                                perturb[1:3]),
        c = card.fit$predict(Cardenolides$b,
                             stats::coef(card.fit$fit),
                             perturb[4:5]),
        l = lma.fit$predict(LMA$b,
                            stats::coef(lma.fit$fit),
                            perturb[6:7])
      )
    )

    detach(obj$vars)
    detach(obj$pars)

    # Initialize with zeros
    H <- array(0, dim=c(n, N))
    C <- array(0, dim=c(n, N))
    L <- array(0, dim=c(n, N))
    cat('done!\n')

    #uses iA indices to store given column of mtx, which has our values for the given variable
    cat('setting herbivory matrix...')
    H[iA] <- mtx[,"w"]
    obj$matrices$H <- H
    cat('done!\n')

    cat('setting cardenolide matrix...')
    C[iA] <- mtx[,"c"]
    obj$matrices$C <- C
    cat('done!\n')

    cat('setting LMA matrix...')
    L[iA] <- mtx[,"l"]
    obj$matrices$L <- L
    cat('done!\n')

    if (update) {
      obj %<>% computeSexualKernel(update = FALSE) %>%
        computeClonalKernel(update = FALSE) %>%
        computeFullKernel()
    }

    if (saveresults) {
      save(obj$matrices$H, obj$matrices$C, obj$matrices$L, file = file.path('data/mwCache',"varDistMats.RData"))
    }

  } else {
    load(file.path('data/mwCache',"varDistMats.RData"))
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

  cat(paste0("Computing seedling recruitment matrix..."))
  R <- t(t(obj$pars$seedling.fit$predict(obj$vars$h_apical$b,
                                         obj$pars$seedling.fit$fit$estimate,
                                         perturb)))%*%t(rep(1,obj$N))

  obj$matrices$R <- R
  cat(paste0("done!\n"))

  if (update) {
    obj %<>% computeSexualKernel()
    #obj <- computeMPM(obj)
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
  attach(obj$vars, warn.conflicts = FALSE)
  attach(obj$pars, warn.conflicts = FALSE)
  attach(obj$matrices, warn.conflicts = FALSE)

  cat(paste0("Compiling sexual kernal..."))
  Ks <- (seedling.emergence + perturb[1])*(seeds.per.pod+perturb[2])*R%*%(P*S*obj$matrices$F)%*%(H*C*L)*herb_avg$dx*h_apical$dx*Cardenolides$dx*LMA$dx
  # image.plot(h_apical, h_apical, t(Ks), col=topo.colors(100))

  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)

  obj$kernels$Ks <- Ks
  cat(paste0("done!\n"))

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

  cat(paste0("Compiling clonal kernal..."))
  if (obj$mdlargs$input == 'meanonly') {          # HERB NEEDS TO BE REPLACED WITH CARDENOLIDES IF WE WANT THIS METHOD TO WORK
    mean.herb_avg <- sum(dist.herb*herb_avg$x)*herb_avg$dx
    mean.buds.per.stem <- budlings.per.stem.fit$predict(mean.herb_avg,
                                                        stats::coef(budlings.per.stem.fit$fit),
                                                        perturb[1:2])
  } else {
    mean.buds.per.stem <- t(budlings.per.stem.fit$predict(Cardenolides$x,
                                                          stats::coef(budlings.per.stem.fit$fit),
                                                          perturb[1:2])) %*%
                          t(t(card.fit$predict(Cardenolides$b,
                                               stats::coef(card.fit$fit),
                                               perturb[1:2])))*Cardenolides$dx
  }

  Kc <- t(t(mean.buds.per.stem*budling.fit$predict(h_apical$b,
                                                   budling.fit$fit$estimate,
                                                   perturb[3:4]))) %*% t(rep(1,obj$N))*h_apical$dx

  detach(obj$vars)
  detach(obj$pars)
  detach(obj$matrices)

  obj$kernels$Kc <- Kc
  cat(paste0("done!\n"))

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

  cat(paste0("Merging kernals..."))
  obj$kernels$K <- obj$kernels$Ks + obj$kernels$Kc
  cat(paste0("done!\n"))
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
        setPodsMatrix(update = FALSE) %>%
        setVarDistMatrices(update=FALSE) %>%
        setSeedlingRecruitmentMatrix(update = FALSE)

      obj %<>%
        computeSexualKernel(update = FALSE) %>%
        computeClonalKernel(update = FALSE) %>%
        computeFullKernel()

      #obj %<>% computeMPM() #commented out by Soren for testing
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

  Kss <- seedling.emergence*seeds.per.pod*(P*G*S*F)%*%H*herb_avg$dx*h_apical.next$dx*h_apical$dx
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
  elas <- matrix(as.vector(sens)*as.vector(K)/lam,nrow=obj$N)
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
  if (!file.exists(file.path('data/mwCache',"parameterAnalysis.RData")) | (compute)) {
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
                                   "Cardenolides")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Survival" %in% params)) {
      cat("Analyzing survival fit...")
      survival_func <- function(x) {obj %>% setSurvivalMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(survival_func, rep(0,5)),
                          pars = obj$pars$surv.fit$pars$unscaled[obj$site,],
                          type = as.character("Survival"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "Cardenolides",
                                   "LMA")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("Pods" %in% params)) {
      cat("Analyzing pods fit...")
      pods_func <- function(x) {obj %>% setPodsMatrix(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(pods_func, rep(0,5)),
                          pars = obj$pars$pods.fit$pars$unscaled[obj$site,],
                          type = as.character("Pods"),
                          name = c("(Intercept)",
                                   "h_apical",
                                   "herb_avg",
                                   "Cardenolides",
                                   "LMA")
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
                          pars = c(obj$pars$seedling.emergence,
                                   obj$pars$seeds.per.pod),
                          type = as.character("Seedlings", "Pods"),
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
                          type = c(rep("Clonal",2), rep("Budlings",2)),
                          name = c("a",
                                   "b",
                                   "mean",
                                   "sd")
        )
        )
      )
      cat("done!\n")
    }

    if (("all" %in% params) | ("VarDist" %in% params)) { ## if there are problems in pertubring just some of the parameters, than check here and verify the whole pipeline of what should be updated and when
      cat("Analyzing variable distributions...")
      varDist_func <- function(x) {obj %>% setVarDistMatrices(perturb = x) %>% analyzeGrowthRate()}
      analysis %<>% bind_rows(
        tbl_df(data.frame(sensitivity = numDeriv::grad(varDist_func, rep(0,7)),
                          pars = c(obj$pars$munched.fit$pmunch,
                                   ParsToMoms(x = obj$pars$munched.fit$fit$estimate,
                                              type = ifelse(!distpars,
                                                            obj$pars$munched.fit$fit$distname,
                                                            "id")),
                                   obj$pars$card.fit$fit$estimate,
                                   obj$pars$lma.fit$fit$estimate),
                          type = c(rep("Herbivory",3), rep("Cardenolides",2), rep("LMA",2)),
                          name = c("pmunch",
                                   "mean",
                                   "sd",
                                   "mean",
                                   "sd",
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
      save(analysis, file = file.path('data/mwCache',"parameterAnalysis.RData"))
    }
  } else {
    load(file.path('data/mwCache',"parameterAnalysis.RData"))
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
ParsToMoms <- function(x, type = "lnorm") {
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
jacParsToMoms <- function(x, type = "lnorm") {
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
MomsToPars <- function(y, type = "lnorm") {
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
jacMomsToPars <- function(y, type = "lnorm") {
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
