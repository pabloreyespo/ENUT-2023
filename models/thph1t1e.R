rm(list = ls())  # .rs.restartR()
source("utils.R")
source("apollo_jaradiaz.R")
set.seed(42)

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

get_data()
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

  P[["model"]] <- apollo_jaradiaz_2pi(jaradiaz_settings = jaradiaz_settings, functionality = functionality)
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


apollo_inputs = apollo_validateInputs(
    apollo_control = model$apollo_control,
    database = database,
    apollo_beta  = model$estimate,
    apollo_fixed = apollo_fixed,
    silent       = T)

pred = apollo_prediction(
  best_model,
  apollo_probabilities,
  apollo_inputs = apollo_inputs,
)

database[,"Tw_pred"] <- pred[,3]

apollo_inputs = apollo_validateInputs(
    apollo_control = model$apollo_control,
    database = database %>% mutate(ec = ec*0.9),
    apollo_beta  = model$estimate,
    apollo_fixed = apollo_fixed,
    silent       = T)

pred = apollo_prediction(
  best_model,
  apollo_probabilities,
  apollo_inputs = apollo_inputs,
)

database[,"Tw_pred"] <- pred[,3]
database %>% # %>% filter((Tw_pred < 168) & (Tw_pred > 0))
  ggplot(aes(x = Tw, y = Tw_pred)) +
  geom_point() +
  geom_smooth()

