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

###########################################################################
###########                     PRE ESTIMACIÓN                   ##########
###########################################################################

nClass <- 2
nvals <- 200
EM <- T
EMiterMax <- 3
penalization <- 1.5
modelName <- paste0("lc_ap_hobbies_", penalization)
socioecon   <- c("female", "mayor45", "underage_in_household", "only_worker", "university", "metropolitana")
especificas <- c("asc", paste0("x_",socioecon))
get_data(
  free_activities = list(
    Tw = c("t_to"),
    Tf1 = c("t_vsyo_csar"),
    Tf2=c("t_vsyo_aa"),
    Tf3=c("t_mcm_leer", "t_mcm_audio", "t_mcm_video", "t_mcm_computador"),
    Tf4 = c("t_cpag_comer")),
  free_expenditures = list(Ef1 = c("alimentos","recreacion","restaurantes","comunicaciones","vestimenta")), baskets = T)
nClass = length(unique(model_data$class))
model_data = model_data %>% filter(ec> 0)
corr_mix <- list()
ntimes <- length(times) - 2
nexpenditures <- max(0, length(expenditures) - 2)
nequations <- ntimes + nexpenditures
for (s in 1:nClass) {corr_mix[[s]] = to_vec(for (i in 1:nequations) for (j in i:nequations) if (i!=j) act_mix[[s]][i]&act_mix[[s]][j] )}

if (ntimes >= 2) { mod_thetas <- to_vec(for (i in times[2:ntimes]) paste0("theta_" ,substring(i,3))) } else {mod_gammas <- c()}
if (nexpenditures >= 1) { mod_phis <- to_vec(for (i in expenditures[1:nexpenditures]) paste0("phi_", substring(i, 3))) } else {mod_deltas <- c()}
mod_sigmas <- to_vec(for (i in 1:nequations) paste0("sigma_",i))
mod_rhos   <- to_vec(for (i in 1:nequations) for (j in i:nequations) if (i!=j) paste0("rho_",i,j))

testvals <- generate_initials_multi_thph(def_sigma = 30, num = nvals , nClass = nClass, especificas = especificas, corners = T,
                                         guess_PH = rep(0.3510, nClass), guess_theta_w = rep(-0.1766, nClass), guess_certainty = 1)
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

apollo_lcPars <- function(apollo_beta, apollo_inputs){
  lcpars <- c()
  lcpars[["PH"]] <- to_list(for (i in 1:nClass) get(paste0("PH_", i)))
  lcpars[["theta_w"]]  <- to_list(for (i in 1:nClass) get(paste0("theta_w_", i)))
  lcpars[["theta_1"]]  <- to_list(for (i in 1:nClass) get(paste0("theta_1_", i)))
  lcpars[["theta_2"]]  <- to_list(for (i in 1:nClass) get(paste0("theta_2_", i)))
  lcpars[["theta_3"]]  <- to_list(for (i in 1:nClass) get(paste0("theta_3_", i)))
  lcpars[["sigma_1"]]  <- to_list(for (i in 1:nClass) get(paste0("sigma_1_", i)))
  lcpars[["sigma_2"]]  <- to_list(for (i in 1:nClass) get(paste0("sigma_2_", i)))
  lcpars[["sigma_3"]]  <- to_list(for (i in 1:nClass) get(paste0("sigma_3_", i)))
  lcpars[["sigma_4"]]  <- to_list(for (i in 1:nClass) get(paste0("sigma_4_", i)))
  lcpars[["rho_12"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_12_", i)))
  lcpars[["rho_13"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_13_", i)))
  lcpars[["rho_14"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_14_", i)))
  lcpars[["rho_23"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_23_", i)))
  lcpars[["rho_24"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_24_", i)))
  lcpars[["rho_34"]]  <- to_list(for (i in 1:nClass) get(paste0("rho_34_", i)))

  V <- list()  #### MOVER A OTRA VARIABLE PARA QUE EL MODELO TENG ALA LIBERTAD DE QUE SEAN NEGATIVAS
  V[["class_1"]] = asc_1 + female * x_female_1 + mayor45 * x_mayor45_1 + underage_in_household * x_underage_in_household_1+ only_worker * x_only_worker_1 + university * x_university_1 +  metropolitana * x_metropolitana_1
  V[["class_2"]] = asc_2 + female * x_female_2 + mayor45 * x_mayor45_2 + underage_in_household * x_underage_in_household_2+ only_worker * x_only_worker_2 + university * x_university_2 +  metropolitana * x_metropolitana_2

  classAlloc_settings <- list(
    classes   = c(class_1 = 1, class_2 = 2),
    utilities = V)

  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()

  for (s in 1:length(pi_values)) {
    mix <- act_mix[[s]]
    obs_times        <- apollo_inputs$database[, times, drop = FALSE]
    obs_expenditures <- apollo_inputs$database[, expenditures, drop = FALSE]
    obs_times_model  <- obs_times[,mix]

    est_sigmas <- to_vec(for (i in 1:nequations) get(paste0("sigma_",i))[[s]])
    est_sigmas <- est_sigmas[mix[1:ntimes]]

    if (ntimes >= 2)        {
      est_thetas <- to_vec(for (th in mod_thetas) get(th)[[s]])
      est_thetas <-  est_thetas[mix[2:ntimes]]
    } else {est_thetas <- c()}
    if (nexpenditures >= 1) { est_phis   <- to_vec(for (ph in mod_phis ) get(ph)[[s]])} else {est_phis <- c()}

    est_corr <- matrix(1, nequations, nequations)
    for (i in 1:(nequations)) for (j in i:(nequations)) if (i != j) { est_corr[j,i] <- est_corr[i,j] <- get(paste0("rho_",i,j))[[s]] }
    est_corr <- est_corr[mix[1:ntimes], mix[1:ntimes]]

    jaradiaz_settings <- list(
      obs_times        = obs_times_model, # number of free time allocations
      obs_expenditures = obs_expenditures, # number of free expednitures
      nequations       = nequations - sum(!mix),
      ntimes           = ntimes - sum(!mix) ,
      nexpenditures    = nexpenditures,
      w                = w,
      ta               = ta,
      TH               = 1,
      PH               = PH[[s]],
      theta_w          = theta_w[[s]],
      thetas           = est_thetas,
      phis             = est_phis,
      sigma            = est_sigmas,
      corr             = est_corr,
      componentName = paste0("class_",s))

    if (sum(!mix)>0) {
      notdoing <- obs_times[, !mix, drop = F]
      pen      <- rowSums(exp(-penalization*(notdoing)))
    } else {
      pen      <- 1
    }

    P[[paste0("class_",s)]]  <- apollo_jaradiaz_2pi(jaradiaz_settings = jaradiaz_settings, functionality = functionality) * pen
  }
  lc_settings <- list(inClassProb = P, classProb=pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
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

  apollo_fixed <- paste0(especificas, "_", 1)
  for (s in 1:nClass) if (sum(!act_mix[[s]]) > 0) {
    apollo_fixed = c(apollo_fixed,
                      paste0(mod_thetas[!act_mix[[s]][2:ntimes]],"_",s),
                      paste0(mod_sigmas[1:ntimes][!act_mix[[s]]],"_",s),
                      paste0(mod_rhos[!corr_mix[[s]]],"_",s))
    apollo_beta[paste0(mod_sigmas[1:ntimes][!act_mix[[s]]],"_",s)] <- 0
  }

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

  model <- NULL
  suppressWarnings({ try({
      if (EM) { #
        invisible(capture.output({
           model <- apollo_lcEM(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                               estimate_settings =  list(
                                 writeIter = FALSE,
                                 silent = T,
                                 scaleHessian = F,
                                 scaleAfterConvergence = F,
                                 validateGrad  = FALSE),
                               lcEM_settings = list(
                                 EMmaxIterations = EMiterMax,
                                 silent =T,
                                 postEM = 0))
        }))
        apollo_beta <- model$estimate
      }
    model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings= est_set)
    })
  })
  try(ans[j,paste0(names(apollo_beta), "_est")]     <-  model$estimate[names(apollo_beta)])
  try(ans[j,paste0(names(apollo_beta), "_estSE")]   <-  model$se[names(apollo_beta)])
  try(ans[j, 'LL']       <-  model$maximum)
  try(ans[j, 'eigen']    <-  model$eigValue)
  try(ans[j, 'state']    <-  model$message)
  try(ans[j, 'code']     <-  model$code)
  try({
    if ((model$code==4&est_set$estimationRoutine=="bgw"|model$code==0&est_set$estimationRoutine=="bfgs") & (model$eigValue < 0)& (model$maximum > best_LL)) { #&
      best_LL <- model$maximum
      best_model <- model
      apollo_saveOutput(best_model)
      tryCatch(apollo_modelOutput(model), error=function(e) NULL)
      #break
    }})
}

name <- paste0("results/", modelName, ".csv")
apollo_saveOutput(best_model)
write.csv(ans, name)

###########################################################################
###########                    POST ESTIMACIÓN                   ##########
###########################################################################









