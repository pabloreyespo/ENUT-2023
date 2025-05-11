rm(list = ls())  # .rs.restartR()
source("utils.R")

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
set.seed(42)

###########################################################################
###########                     PRE ESTIMACIÓN                   ##########
###########################################################################
nClass <- 4
nvals <- 200
EM <- T
EMiterMax <- 3
modelName <- "4C-quadratic"
socioecon   <- c("female", "mayor45", "underage_in_household", "only_worker", "university", "metropolitana")
especificas <- c("asc", paste0("x_",socioecon))
get_data(
  free_activities = list(
    Tw = c("t_to"),
    Tfleisure = c('t_vsyo_csar', 't_vsyo_aa', 't_mcm_leer', 't_mcm_audio', 't_mcm_video', 't_mcm_computador')),
  free_expenditures = list(Ef1 = c("alimentos","recreacion","restaurantes","comunicaciones","vestimenta")))

model_data = model_data %>% mutate(wbar = mean(w), ta= 168 - Tc, tabar = mean(168-Tc), Ecbar = mean(Ec))
testvals <- generate_initials_quadratic(def_sigma = 100, num = nvals, nClass = nClass, especificas = especificas, guess_certainty = 0.01)
ans <- c(paste0(names(testvals), '_initial'),
         paste0(names(testvals), '_est'),
         paste0(names(testvals), '_estSE'),
         'LL',"eigen","state","code")

ans <- matrix(NA, nrow=nvals, ncol=length(ans), dimnames=list(NULL, ans))
k <- length(testvals)

################################################
###########                     ESTIMACION                       ##########
###########################################################################
apollo_initialise()
database <- model_data

apollo_lcPars <- function(apollo_beta, apollo_inputs){
  lcpars <- c()
  lcpars[["alpha"]] = to_list(for (i in 1:nClass) get(paste0("alpha_", i)))
  lcpars[["beta_w1"]] = to_list(for (i in 1:nClass) get(paste0("beta_w1_", i)))
  lcpars[["beta_ta1"]] = to_list(for (i in 1:nClass) get(paste0("beta_ta1_", i)))
  lcpars[["beta_Ec1"]] = to_list(for (i in 1:nClass) get(paste0("beta_Ec1_", i)))
  lcpars[["beta_w2"]] = to_list(for (i in 1:nClass) get(paste0("beta_w2_", i)))
  lcpars[["beta_ta2"]] = to_list(for (i in 1:nClass) get(paste0("beta_ta2_", i)))
  lcpars[["beta_Ec2"]] = to_list(for (i in 1:nClass) get(paste0("beta_Ec2_", i)))
  lcpars[["beta_wta"]] = to_list(for (i in 1:nClass) get(paste0("beta_wta_", i)))
  lcpars[["beta_wEc"]] = to_list(for (i in 1:nClass) get(paste0("beta_wEc_", i)))
  lcpars[["beta_taEc"]] = to_list(for (i in 1:nClass) get(paste0("beta_taEc_", i)))
  lcpars[["sigma"]] <- to_list(for (i in 1:nClass) get(paste0("sigma_", i)))
  V <- list()  #### MOVER A OTRA VARIABLE PARA QUE EL MODELO TENG ALA LIBERTAD DE QUE SEAN NEGATIVAS
  V[["class_1"]] = asc_1 + female * x_female_1 + mayor45 * x_mayor45_1 + underage_in_household * x_underage_in_household_1+ only_worker * x_only_worker_1 + university * x_university_1 +  metropolitana * x_metropolitana_1
  V[["class_2"]] = asc_2 + female * x_female_2 + mayor45 * x_mayor45_2 + underage_in_household * x_underage_in_household_2+ only_worker * x_only_worker_2 + university * x_university_2 +  metropolitana * x_metropolitana_2
  V[["class_3"]] = asc_3 + female * x_female_3 + mayor45 * x_mayor45_3 + underage_in_household * x_underage_in_household_3+ only_worker * x_only_worker_3 + university * x_university_3 +  metropolitana * x_metropolitana_3
  V[["class_4"]] = asc_4 + female * x_female_4 + mayor45 * x_mayor45_4 + underage_in_household * x_underage_in_household_4+ only_worker * x_only_worker_4 + university * x_university_4 +  metropolitana * x_metropolitana_4

  classAlloc_settings <- list(
    classes   = c(class_1 = 1, class_2 = 2, class_3 = 3, class_4 = 4),
    utilities = V)
  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()

  for (s in 1:length(pi_values)) {
    normalDensity_settings <- list(
     outcomeNormal   = Tw,
     xNormal         = alpha[[s]] + beta_w1[[s]]*(w-wbar) + beta_w2[[s]]*(w-wbar)^2 + beta_ta1[[s]]*(ta-tabar) + beta_ta2[[s]]*(ta-tabar)^2 +  beta_Ec1[[s]]*(Ec-Ecbar) + beta_Ec2[[s]]*(Ec-Ecbar)^2 +
                        beta_wta[[s]]*(w-wbar)*(ta-tabar) + beta_wEc[[s]]*(w-wbar)*(Ec-Ecbar) + beta_taEc[[s]]*(ta-tabar)*(Ec-Ecbar),
     mu    = 0, sigma = sigma[[s]],
     componentName = paste0("class_",s))

    P[[paste0("class_",s)]]  <- apollo_normalDensity(normalDensity_settings, functionality = functionality) }
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
      modelDescr      = "quadratic labor supply",
      indivID         = "id_persona",
      outputDirectory = "output")
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

out <- post_eval_latent_class(modelName, model_data, set_class = T)
best_model <- out$model
model_data <- out$data
apollo_modelOutput(best_model)

cat("probabilistic assign:", "\n")
cat(colSums(model_data[, paste0("pi_",1:nClass)]) / nrow(model_data), "\n")
cat("deterministic assign:", "\n")
cat(model_data[,paste0("pi_",1:nClass)] %>% max.col() %>% table() / nrow(model_data), "\n")

pi <- model_data %>% group_by(class) %>%
  summarise(pi_1_prom = mean(pi_1),
            pi_2_prom = mean(pi_2),
            pi_3_prom = mean(pi_3),
            pi_4_prom = mean(pi_4))

print(pi)