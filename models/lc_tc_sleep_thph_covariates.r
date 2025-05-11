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
modelName <- "Tc-sleep-ENUT-THPH1T1E-covariates"
socioecon   <- c("female", "mayor45", "underage_in_household", "only_worker", "university", "metropolitana")
especificas <- c("asc", paste0("x_",socioecon))
get_data_tc(especificas = socioecon, disputed = "t_sleep")
model_data = model_data %>% filter(ec_1> 0, ec_2> 0) # , ec_1 < 1, ec_2 < 1
covariates <- c("female", "mayor45", "underage_in_household", "only_worker")
testvals <- generate_initials_simple_tc_thph(def_sigma = 20, num = nvals, especificas = especificas,
                                    guess_PH = c(0.7180,0.7180), guess_theta_w = c(-0.1578, -0.1578) , guess_certainty = 0.5,
                                             covariates = covariates)
ans <- c(paste0(names(testvals), '_initial'),
         paste0(names(testvals), '_est'),
         paste0(names(testvals), '_estSE'),
         'LL',"eigen","state","code")
ans <- matrix(NA, nrow=nvals, ncol=length(ans), dimnames=list(NULL, ans))
k <- length(testvals)
###########################################################################
###########                     ESTIMACION                       ##########
###########################################################################
apollo_initialise()
database <- model_data
apollo_lcPars <- function(apollo_beta, apollo_inputs){
  lcpars <- c()
  lcpars[["PH"]] <- to_list(for (i in 1:nClass) get(paste0("PH_", i)))
  for (co in covariates) { lcpars[[paste0("PH_",co)]] <- to_list(for (i in 1:nClass) get(paste0("PH_",co,"_", i)))}
  lcpars[["theta_w"]]  <- to_list(for (i in 1:nClass) get(paste0("theta_w_", i)))
  for (co in covariates) { lcpars[[paste0("theta_w_",co)]] <- to_list(for (i in 1:nClass) get(paste0("theta_w_",co,"_", i)))}
  lcpars[["sigma"]] <- to_list(for (i in 1:nClass) get(paste0("sigma_", i)))
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
    PH_value       = PH[[s]]      +female*PH_female[[s]]      +mayor45*PH_mayor45[[s]]       + underage_in_household * PH_underage_in_household[[s]] + only_worker * PH_only_worker[[s]]              #+ university * PH_university[[s]]      #+ metropolitana * PH_metropolitana[[s]]                       #
    theta_w_value  = theta_w[[s]] +female*theta_w_female[[s]] +mayor45*theta_w_mayor45[[s]]  + underage_in_household * theta_w_underage_in_household[[s]] + only_worker  * theta_w_only_worker[[s]]   #+ university * theta_w_university[[s]] #+ metropolitana * theta_w_metropolitana[[s]] #
    jaradiaz_settings <- list(Tw          = apollo_inputs$database[, paste0("Tw_",s)],
                              Tc          = apollo_inputs$database[, paste0("Tc_",s)],
                              Ec          = apollo_inputs$database[, paste0("Ec_",s)],
                              w           = w,
                              ta          = ta,
                              TH          = 1,
                              PH          = PH_value,
                              theta_w     = theta_w_value,
                              sigma       = sigma[[s]],
                              componentName = paste0("class_",s))
    P[[paste0("class_",s)]] <- apollo_jaradiaz_2pi(jaradiaz_settings = jaradiaz_settings, functionality = functionality)
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
  apollo_beta <- c()
  for (b in names(testvals)) {apollo_beta[b] <- as.numeric(testvals[[b]][j])}

  ans[j, 1:k]  <- apollo_beta
  apollo_fixed <- paste0(especificas, "_", 1)

  apollo_control = list(
    modelName       = modelName,
    modelDescr      = "Jara-Diaz model",
    indivID         = "id_persona",
    outputDirectory = "output",
    analyticGrad = FALSE,
    noValidation = FALSE,
    noDiagnostics = TRUE,
    workInLogs    = FALSE)
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
  try(ans[j, 'LL']                                  <-  model$maximum)
  try(ans[j, 'eigen']                               <-  model$eigValue)
  try(ans[j, 'state']                               <-  model$message)
  try(ans[j, 'code']                                <-  model$code)
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
