rm(list = ls())  # .rs.restartR()
pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats", "lubridate", "haven", "reshape2","apollo", "comprehenr", "Hmisc")
est_set <- list(
  writeIter = FALSE,
  silent = T,
  maxIterations=500,
  scaleHessian = F,
  scaleAfterConvergence = F,
  estimationRoutine = "bgw",
  hessianRoutine = "maxLik",
  bgw_settings = list(maxFunctionEvals = 1000),
  validateGrad  = FALSE
)
invisible(lapply(pkgs, library, character.only=TRUE))
source("utils.R")
source("apollo_jaradiaz.R")
set.seed(42)

###############################################################################################################

##############################################################################################################


get_data()
get_data(
  free_activities = list(
    Tw = c("t_to"),
    Tfleisure = c('t_vsyo_csar', 't_vsyo_aa', 't_mcm_leer', 't_mcm_audio', 't_mcm_video', 't_mcm_computador')),
  free_expenditures = list(Ef1 = c("alimentos","recreacion","restaurantes","comunicaciones","vestimenta")))

model_data <- model_data %>% filter(ec>0)
nvals <- 100
modelName <- "ENUT-THPH1T1E"
testvals <- generate_initials_simple_thph(def_sigma = 10, num = nvals, nClass = 1)
ans <- c(paste0(names(testvals), '_initial'),
         paste0(names(testvals), '_est'),
         paste0(names(testvals), '_estSE'),
         'LL','eigen',"state","code")

ans <- matrix(NA, nrow=nvals, ncol=length(ans), dimnames=list(NULL, ans))
k <- length(testvals)
###########################################################################
###########                     ESTIMACION                       ##########
###########################################################################
apollo_initialise()
database <- model_data
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()

  jaradiaz_settings <- list(
    Tw          = apollo_inputs$database[, "Tw"],
    Tc          = apollo_inputs$database[, "Tc"],
    Ec          = apollo_inputs$database[, "Ec"],
    w           = w,
    ta          = ta,
    TH          = 1,
    PH          = PH,
    theta_w     = theta_w,
    sigma       = sigma,
    componentName = "model")

  P[["model"]] <- apollo_jaradiaz(jaradiaz_settings = jaradiaz_settings, functionality = functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

best_model <- NULL
best_LL    <- -Inf

for (j in 1:nvals) {
  cat("EVALUANDO",j,"/",nvals,"...\n")

  apollo_beta = c()
  for (b in names(testvals)) {apollo_beta[b] = as.numeric(testvals[[b]][j])}
  ans[j, 1:k] = apollo_beta

  apollo_fixed <- c()

  apollo_inputs = apollo_validateInputs(
    apollo_control = list(
      modelName       = modelName,
      modelDescr      = "Modelo Jara-Diaz",
      indivID         = "id_persona",
      outputDirectory = "output",
      noValidation = TRUE,
      analyticGrad = TRUE,
      noDiagnostics = TRUE,
      workInLogs    = FALSE,
      debug         = F),
    apollo_beta  = apollo_beta,
    apollo_fixed = apollo_fixed,
    silent       = T)

  model <- NULL
  suppressWarnings({
    try(model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities,
                          apollo_inputs, estimate_settings=est_set))})

  try(ans[j,paste0(names(testvals), '_est')]   <- model$estimate)
  try(ans[j,paste0(names(testvals), '_estSE')] <- model$se)
  try(ans[j, 'LL']       <-  model$maximum)
  try(ans[j, 'eigen']    <-  model$eigValue)
  try(ans[j, 'state']    <-  model$message)
  try(ans[j, 'code']     <-  model$code)
  try({
    if ((model$code==4&est_set$estimationRoutine=="bgw"|model$code==0&est_set$estimationRoutine=="bfgs") & (model$eigValue < 0)& (model$maximum > best_LL)) { #&
      best_LL <- model$maximum
      best_model <- model
      tryCatch(apollo_modelOutput(model), error=function(e) NULL)
      break
    }})
}

library(glue)
report_values_of_time_thph(best_model, model_data)




