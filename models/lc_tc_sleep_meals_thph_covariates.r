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
nClass <- 4
nvals <- 500
EM <- T
EMiterMax <- 3
modelName <- "Tc-sleep-meals-ENUT-THPH1T1E"
socioecon   <- c("female", "mayor45", "underage_in_household", "only_worker", "university", "metropolitana")
especificas <- c("asc", paste0("x_",socioecon))
get_data_2tc(especificas = socioecon, disputed = c("t_sleep", "t_meals"))
model_data = model_data %>% filter(ec_1> 0, ec_2> 0, ec_3> 0, ec_4> 0)
covariates <- c("female", "mayor45", "underage_in_household", "only_worker")
guess_PH = rep(0.7180, nClass)
guess_theta_w = rep(-0.1578, nClass)
testvals <- generate_initials_simple_tc_thph(
  def_sigma = 20, 
  num = nvals, 
  especificas = especificas, 
  nClass=nClass,
  guess_PH = guess_PH, 
  guess_theta_w = guess_theta_w, 
  guess_certainty = 0.3,
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
  V[["class_1"]]= asc_1+ female*x_female_1+ mayor45*x_mayor45_1+ underage_in_household*x_underage_in_household_1+ only_worker*x_only_worker_1+ university*x_university_1+ metropolitana*x_metropolitana_1
  V[["class_2"]]= asc_2+ female*x_female_2+ mayor45*x_mayor45_2+ underage_in_household*x_underage_in_household_2+ only_worker*x_only_worker_2+ university*x_university_2+ metropolitana*x_metropolitana_2
  V[["class_3"]]= asc_3+ female*x_female_3+ mayor45*x_mayor45_3+ underage_in_household*x_underage_in_household_3+ only_worker*x_only_worker_3+ university*x_university_3+ metropolitana*x_metropolitana_3
  V[["class_4"]]= asc_4+ female*x_female_4+ mayor45*x_mayor45_4+ underage_in_household*x_underage_in_household_4+ only_worker*x_only_worker_4+ university*x_university_4+ metropolitana*x_metropolitana_4

  classAlloc_settings <- list(
    classes   = c(class_1= 1, class_2= 2, class_3= 3, class_4= 4),
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

###########################################################################
###########                    POST ESTIMACIÓN                   ##########
###########################################################################
name <- paste0("results/", modelName, ".csv")
best_model <- apollo_loadModel(paste0("output/", modelName))
ans_out <- read.csv(name)
apollo_modelOutput(best_model)
get_membership_probabilities(best_model)
report_values_of_time_classes_thph(best_model, dbs = dbs, tc = T)
cat("probabilistic assign:", "\n")
cat(colSums(dbs[, paste0("pi_",1:nClass)]) / nrow(dbs), "\n")
cat("deterministic assign:", "\n")
cat(dbs[,paste0("pi_",1:nClass)] %>% max.col() %>% table() / nrow(dbs), "\n")
proportions <- c(socioecon,"menor25","n_menores6","n_menores12","vive_pareja","e0", "e1", "e2", "e3","q1", "q2", "q3", "q4", "q5",
                 "Norte_grande", "Norte_chico", "Zona_centro", "Zona_sur", "Zona_austral",
                 "fuentes_externas", "trabajador_obrero", "trabajador_privilegiado")
socioecon_plus   <- c("edad_anios", "n_personas", "n_menores", "w", "I","ec_1","ec_2","ing_trab")
mean_values <- get_segment_probabilities(dbs, socioecon_plus, proportions, tc=T)
print(mean_values)

plot_thph_model(best_model, dbs ,tc=T)
plot_4tc_difference_covariates(models, lc_model, str_to_title(act))
plot_4tc_difference_covariates_vot(models, lc_model, model_data, str_to_title(act))