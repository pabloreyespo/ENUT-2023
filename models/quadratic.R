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
source('utils.R')
set.seed(42)

###############################################################################################################

##############################################################################################################

get_data()
model_data = model_data %>% mutate(wbar = mean(w), ta= 168 - Tc, tabar = mean(168-Tc), Ecbar = mean(Ec))
#model_data <- type.convert(model_data, as.is =TRUE)

nvals <- 100
modelName <- "ENUT-quadratic"
testvals <- generate_initials_quadratic(def_sigma=50, num = 1, nClass = 1, especificas = NULL, guess_certainty = 0 )
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

  normalDensity_settings <- list(
   outcomeNormal   = Tw,
   xNormal         = alpha + beta_w1*(w-wbar) + beta_w2*(w-wbar)^2 + beta_ta1*(ta-tabar) + beta_ta2*(ta-tabar)^2 + beta_Ec1*(Ec-Ecbar) + beta_Ec2*(Ec-Ecbar)^2 +
                      beta_wta*(w-wbar)*(ta-tabar) + beta_wEc*(w-wbar)*(Ec-Ecbar) + beta_taEc*(ta-tabar)*(Ec-Ecbar),
   mu    = 0,
   sigma = sigma
  )

  P[["model"]] <- apollo_normalDensity(normalDensity_settings, functionality = functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

best_model <- NULL
best_LL    <- -Inf

apollo_beta = c()
for (b in names(testvals)) {apollo_beta[b] = testvals[[b]][1]}
apollo_fixed <- c()

apollo_inputs = apollo_validateInputs(
  apollo_control = list(
    modelName       = modelName,
    modelDescr      = "Modelo Jara-Diaz",
    indivID         = "id_persona",
    outputDirectory = "output"),
  apollo_beta  = apollo_beta,
  apollo_fixed = apollo_fixed,
  silent       = T)

model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities,
                          apollo_inputs, estimate_settings=est_set)
apollo_modelOutput(model)

