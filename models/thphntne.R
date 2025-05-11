rm(list = ls())  # .rs.restartR()
source("utils.R")
source("apollo_jaradiaz.R")

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

diagonal_correlation = F

set.seed(42)

###########################################################################
###########                     PRE ESTIMACIÓN                   ##########
###########################################################################

nvals <- 100
modelName <- "ENUT-THPHNTNE"
get_data(free_activities = list(
            t_paid_work             = c("t_to"),
            t_leisure_socialization = c("t_vsyo_csar"),
            t_leisure_hobbies       = c("t_vsyo_aa") ,
            t_leisure_reading       = c("t_mcm_leer"),
            t_leisure_audio         = c("t_mcm_audio"),
            t_leisure_video         = c("t_mcm_video"),
            t_leisure_computer         = c("t_mcm_computador"),
            t_meals = c("t_cpag_comer")
            ))
model_data = model_data %>% filter(ec>0)

ntimes <- length(times) - 2
nexpenditures <- max(0, length(expenditures) - 2)
nequations <- ntimes + nexpenditures

if (ntimes >= 2) { mod_thetas <- to_vec(for (i in times[2:ntimes]) paste0("theta_" ,substring(i,3))) } else {mod_gammas <- c()}
if (nexpenditures >= 1) { mod_phis <- to_vec(for (i in expenditures[1:nexpenditures]) paste0("phi_", substring(i, 3))) } else {mod_deltas <- c()}
mod_sigmas <- to_vec(for (i in 1:nequations) paste0("sigma_",i))
mod_rhos   <- to_vec(for (i in 1:nequations) for (j in i:nequations) if (i!=j) paste0("rho_",i,j))

testvals <- generate_initials_multi_thph(def_sigma = 30, num = nvals, guess_PH = c(0.51006), guess_theta_w = c(-0.01929))
ans <- c(paste0(names(testvals), '_initial'),
         paste0(names(testvals), '_est'),
         paste0(names(testvals), '_estSE'),
         'LL','eigen',"state","code")

ans <- matrix(NA, nrow=nvals, ncol=length(ans),
              dimnames=list(NULL, ans))

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

  obs_times <- apollo_inputs$database[, times, drop = FALSE]
  obs_expenditures <- apollo_inputs$database[, expenditures, drop = FALSE]

  est_sigmas <- to_vec(for (i in 1:nequations) get(paste0("sigma_",i)))
  if (ntimes >= 2)        { est_thetas <- to_vec(for (th in mod_thetas) get(th))} else {est_thetas <- c()}
  if (nexpenditures >= 1) { est_phis   <- to_vec(for (ph in mod_phis  ) get(ph))} else {est_phis <- c()}

  est_corr <- matrix(1, nequations, nequations)
  for (i in 1:(nequations)) for (j in i:(nequations)) if (i != j) { est_corr[j,i] <- est_corr[i,j] <- get(paste0("rho_",i,j)) }

  jaradiaz_settings <- list(
    obs_times        = obs_times, # number of free time allocations
    obs_expenditures = obs_expenditures, # number of free expednitures
    nequations       = nequations,
    ntimes           = ntimes,
    nexpenditures    = nexpenditures,
    w                = w,
    ta               = ta,
    TH               = 1,
    PH               = PH,
    theta_w          = theta_w,
    thetas           = est_thetas,
    phis             = est_phis,
    sigma            = est_sigmas,
    corr             = est_corr,
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

  apollo_fixed <- c() # mod_rhos
  if (diagonal_correlation) {apollo_fixed = c(mod_rhos)} else {apollo_fixed <- c()}

  apollo_control = list(
      modelName       = modelName,
      modelDescr      = "Jara-Diaz model",
      indivID         = "id_persona",
      outputDirectory = "output",
      analyticGrad = FALSE,
      noValidation = FALSE,
      noDiagnostics = TRUE,
      workInLogs    = FALSE,
      debug         = F)

  apollo_inputs = apollo_validateInputs(
    apollo_control = apollo_control,
    apollo_beta  = apollo_beta,
    apollo_fixed = apollo_fixed,
    silent       = T)

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
    if ((model$code==4&est_set$estimationRoutine=="bgw"|model$code==0&est_set$estimationRoutine=="bfgs") & (model$eigValue <= 0)& (model$maximum > best_LL)) { #&
      best_LL <- model$maximum
      best_model <- model
      tryCatch(apollo_modelOutput(model), error=function(e) NULL)
      break
    }})
}

name <- paste0("results/", modelName, ".csv")
apollo_saveOutput(best_model)

###########################################################################
###########                    POST ESTIMACIÓN                   ##########
###########################################################################

best_model <- apollo_loadModel(paste0("output/", modelName))
apollo_modelOutput(best_model)
report_values_of_time_thph(best_model, model_data)







