getPackages <- function() {
  pkgs <- c("apollo", "comprehenr", "dplyr", "fitdistrplus", "forcats","geomtextpath", "ggplot2", "ggpubr",
            "glue", "gridExtra", "haven", "Hmisc", "lubridate", "matrixcalc", "purrr", "randcorr",
            "RColorBrewer", "readr", "reshape2", "stringr", "tibble", "tidyr", "see")
  invisible(lapply(pkgs, library, character.only=TRUE))
}

colores <- c( "#2D4263", "#005F73","#0A9396", "#D0CA6B","#E9D8A6", "#EE9B00", "#CA6702","#C84B31") #, " tercero # quinto ,
colores <- c( "#2D4263", "#005F73", "#8E9F4C", "#EE9B00","#C84B31") #, " tercero # quinto ,
colores <- c("#0071b2","#009e74","#e69f00","#d55e00")
getPalette <- colorRampPalette(colores)
getPackages()

theme_tesis <- theme(
  strip.text = element_text(color = "#2D4263", size = 12),
  strip.background = element_rect(fill = "white"),
  axis.text = element_text(colour="#2D4263", size = 12),
  legend.title = element_text(colour="#2D4263", size = 13),
  legend.text  = element_text(colour="#2D4263", size = 12),
  axis.title = element_text(colour="#2D4263", size = 12),
  plot.title = element_text(colour="#2D4263", size = 15),
  plot.caption=element_text(colour="#2D4263"),
  legend.position = "bottom")

all_activities   <- c(
  "t_to",
  "t_tcnr_ce",
  "t_tcnr_re",
  "t_tcnr_oac",
  "t_tdnr_psc",
  "t_tdnr_lv",
  "t_tdnr_lrc",
  "t_tdnr_mrm",
  "t_tdnr_admnhog",
  "t_tdnr_comphog",
  "t_tdnr_cmp",
  "t_tvaoh_tv",
  "t_tvaoh_oh",
  "t_cpaf_cp",
  "t_cpag_comer",
  "t_cpag_dormir",
  "t_ed",
  "t_vsyo_csar",
  "t_vsyo_aa",
  "t_mcm_leer",
  "t_mcm_audio",
  "t_mcm_video",
  "t_mcm_computador",
  "t_tt1",
  "t_tt2"
)
all_expenditures <- c(
  "alimentos",
  "recreacion",
  "restaurantes",
  "comunicaciones",
  "vestimenta",
  "cuentas",
  "hogar",
  "salud",
  "transporte",
  "educacion"
)

########################################################################################################################
###################################               LOAD DATABASES               #########################################
########################################################################################################################

get_data <- function(
  free_activities = list(Tw = c("t_to"),
                         Tfleisure = c("t_vsyo_csar", "t_vsyo_aa","t_mcm_leer", "t_mcm_audio", "t_mcm_video", "t_mcm_computador"),
                         Tfmeals = c("t_cpag_comer")),
  free_expenditures = list(Ef1 = c("alimentos","recreacion","restaurantes","comunicaciones"),
                           Ef2 = c("vestimenta")),
  baskets = FALSE
) {

  times <<- c(names(free_activities), 'Tc')
  expenditures <<- c(names(free_expenditures), "Ec")
  fixed_income <- c("ing_jub_aps", "ing_gpp")

  model_data <- haven::read_dta("data/enut-ii-25G.dta") %>%
    filter(edad_años >= 18, es_trabajador == 1, ing_personal > 0, w > 0, total_expenses/ing_personal <= 5) %>%
    mutate(female = sexo,
           menor25 = edad_años <  25,
           mayor45 = edad_años >= 45,
           mayor66 = edad_años >= 60,
           university = nivel_escolaridad ==  'universitaria',
           underage_in_household = case_when(n_menores >= 1 ~ 1, T ~ 0),
           children_in_household = case_when(n_menores6 + n_menores12 >= 1 ~ 1, T ~ 0),
           only_worker = case_when(n_trabajadores == 1 & trabaja == 1 ~ 1, T ~ 0),
           metropolitana = macrozona == "metropolitana",
           norte = macrozona == "norte",
           centro = macrozona == "centro",
           sur = macrozona == "sur"
    )

  model_data[, "ta"] <- 168
  used_activities <- c()
  for (act in names(free_activities)) {
    model_data[, act] <- rowSums(model_data[, free_activities[[act]], drop = F])
    used_activities <- c(used_activities, free_activities[[act]]) }
  model_data[, "Tc"] <- rowSums(model_data[, all_activities[!all_activities %in% used_activities]])

  used_expenditures <- c()
  for (exp in names(free_expenditures)) {
    model_data[, exp] <- rowSums(model_data[, free_expenditures[[exp]], drop = F])
    used_expenditures <- c(used_expenditures, free_expenditures[[exp]]) }
  model_data[, "Ec"] <- rowSums(model_data[, all_expenditures[!all_expenditures %in% used_expenditures]]) - rowSums(model_data[, fixed_income])
  model_data[, "ec"] <- model_data[, "Ec"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc"]))
  model_data["I"] <- rowSums(model_data[, fixed_income])

  times_labeled <- paste0(times, "_label")
  model_data[times_labeled] <- model_data[, times]
  model_data[times_labeled] <- 0
  model_data[,times_labeled][model_data[,times] > 0] <- 1
  model_data[,"class"] <- apply( model_data[ , c(times_labeled) ] , 1 , paste , collapse = "" )

  if (baskets) {
    nClass <- 0
    act_mix <- list()
    classes <- table(model_data$class)
    classes <- names(classes[classes > 200]) #
    model_data = model_data %>% filter(class %in% classes)

    for (cla in unique(model_data$class) %>% sort()){
      nClass <- nClass + 1
      mask <- model_data["class"] == cla
      model_data[mask, "alt"] <- nClass
      act_mix[[nClass]] <- mix <- as.logical(as.integer(strsplit(cla, "")[[1]]))
    }

    nClass <<- nClass
    act_mix <<- act_mix
  }

  model_data <- type.convert(model_data, as.is =TRUE)
  model_data <<- model_data# %>% dplyr::select(id_persona, ing_personal, w, ta, all_of(times), all_of(expenditures), all_of(especificas))
}

get_data_tc <- function(especificas = c(), disputed = "t_sleep") {
  times <<- c('Tw', 'Tfleisure', 'Tc')
  expenditures <<- c("Ef","Ec")
  exp_committed <- c("cuentas", "hogar","salud", "transporte" , "educacion")
  exp_1 <- c("alimentos","recreacion","restaurantes","vestimenta","comunicaciones")
  fixed_income <- c("ing_jub_aps", "ing_gpp")

  model_data <- haven::read_dta("data/enut-ii-11G.dta")  %>%
    filter(edad_años >= 18, es_trabajador == 1, ing_personal > 0, w > 0, total_expenses/ing_personal <= 5) %>%
    mutate(female = sexo,
           menor25 = edad_años <  25,
           mayor45 = edad_años >= 45,
           mayor66 = edad_años >= 60,
           university = nivel_escolaridad ==  'universitaria',
           underage_in_household = case_when(n_menores >= 1 ~ 1, T ~ 0),
           children_in_household = case_when(n_menores6 + n_menores12 >= 1 ~ 1, T ~ 0),
           only_worker = case_when(n_trabajadores == 1 & trabaja == 1 ~ 1, T ~ 0),
           metropolitana = macrozona == "metropolitana",
           norte = macrozona == "norte",
           centro = macrozona == "centro",
           sur = macrozona == "sur"
    )

  leisure_activities <- c("t_leisure", disputed)
  acts  <-c("t_personal_care", "t_sleep","t_meals", "t_commute1", "t_commute2", "t_care_work", "t_domestic_work", "t_unpaid_voluntary", "t_education")
  non_leisure_activities <- acts[!(acts %in% leisure_activities)]
  model_data[, "ta"] <- 168
  model_data[, "Tw_1"]        <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_1"] <- model_data[, "t_leisure"]
  model_data[, "Tc_1"] <- rowSums(model_data[, acts])
  model_data[, "Ec_1"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_1"] <- model_data[, "Ec_1"] / (model_data[, "w"] * (model_data[, "ta"] - model_data[, "Tc_1"]))

  model_data[, "Tw_2"] <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_2"] <- rowSums(model_data[, leisure_activities])
  model_data[, "Tc_2"] <- rowSums(model_data[, non_leisure_activities])
  model_data[, "Ec_2"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_2"] <- model_data[, "Ec_2"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc_2"]))

  model_data[, "Ef_1"] <- rowSums(model_data[, exp_1])
  model_data[, "Ef_2"] <- rowSums(model_data[, exp_1])
  model_data["I"] <- rowSums(model_data[, fixed_income])
  model_data <- type.convert(model_data, as.is =TRUE)
  model_data <<- model_data# %>% dplyr::select(id_persona, ing_personal, w, ta, all_of(paste0(times,"_",1)), all_of(paste0(times,"_",2)), all_of(paste0(expenditures,"_",1)), all_of(paste0(expenditures,"_",2)), all_of(especificas))
}

get_data_2tc <- function(especificas = c(), disputed = c("t_sleep", "t_meals")) {
  times <<- c('Tw', 'Tfleisure', 'Tc')
  expenditures <<- c("Ef","Ec")
  exp_committed <- c("cuentas", "hogar","salud", "transporte" , "educacion")
  exp_1 <- c("alimentos","recreacion","restaurantes","vestimenta","comunicaciones")
  fixed_income <- c("ing_jub_aps", "ing_gpp")

  model_data <- haven::read_dta("data/enut-ii-11G.dta")  %>%
    filter(edad_años >= 18, es_trabajador == 1, ing_personal > 0, w > 0, total_expenses/ing_personal <= 5) %>%
    mutate(female = sexo,
           menor25 = edad_años <  25,
           mayor45 = edad_años >= 45,
           mayor66 = edad_años >= 60,
           university = nivel_escolaridad ==  'universitaria',
           underage_in_household = case_when(n_menores >= 1 ~ 1, T ~ 0),
           children_in_household = case_when(n_menores6 + n_menores12 >= 1 ~ 1, T ~ 0),
           only_worker = case_when(n_trabajadores == 1 & trabaja == 1 ~ 1, T ~ 0),
           metropolitana = macrozona == "metropolitana",
           norte = macrozona == "norte",
           centro = macrozona == "centro",
           sur = macrozona == "sur"
    )

  leisure_activities <- c("t_leisure", disputed)
  acts  <-c("t_personal_care", "t_sleep","t_meals", "t_commute1", "t_commute2", "t_care_work", "t_domestic_work", "t_unpaid_voluntary", "t_education")
  non_leisure_activities <- acts[!(acts %in% leisure_activities)]
  model_data[, "ta"] <- 168

  model_data[, "Tw_1"]        <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_1"] <- model_data[, "t_leisure"]
  model_data[, "Tc_1"] <- rowSums(model_data[, acts])
  model_data[, "Ec_1"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_1"] <- model_data[, "Ec_1"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc_1"]))

  model_data[, "Tw_2"] <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_2"] <- rowSums(model_data[, leisure_activities]) - model_data[, disputed[1]]
  model_data[, "Tc_2"] <- rowSums(model_data[, non_leisure_activities]) + model_data[, disputed[1]]
  model_data[, "Ec_2"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_2"] <- model_data[, "Ec_2"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc_2"]))

  model_data[, "Tw_3"] <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_3"] <- rowSums(model_data[, leisure_activities]) - model_data[, disputed[2]]
  model_data[, "Tc_3"] <- rowSums(model_data[, non_leisure_activities]) + model_data[, disputed[2]]
  model_data[, "Ec_3"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_3"] <- model_data[, "Ec_3"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc_3"]))

  model_data[, "Tw_4"] <- model_data[, "t_paid_work"]
  model_data[, "Tfleisure_4"] <- rowSums(model_data[, leisure_activities])
  model_data[, "Tc_4"] <- rowSums(model_data[, non_leisure_activities])
  model_data[, "Ec_4"] <- rowSums(model_data[, exp_committed]) - rowSums(model_data[, fixed_income])
  model_data[, "ec_4"] <- model_data[, "Ec_4"] / (model_data[, "w"] * (model_data[, "ta"]-model_data[, "Tc_4"]))

  model_data[, "Ef_1"] <- rowSums(model_data[, exp_1])
  model_data[, "Ef_2"] <- rowSums(model_data[, exp_1])
  model_data[, "Ef_3"] <- rowSums(model_data[, exp_1])
  model_data[, "Ef_4"] <- rowSums(model_data[, exp_1])
  model_data["I"] <- rowSums(model_data[, fixed_income])
  model_data <- type.convert(model_data, as.is =TRUE)
  model_data <<- model_data# %>% dplyr::select(id_persona, ing_personal, w, ta, all_of(paste0(times,"_",1)), all_of(paste0(times,"_",2)), all_of(paste0(expenditures,"_",1)), all_of(paste0(expenditures,"_",2)), all_of(especificas))
}

########################################################################################################################
###################################               INTIIAL VALUES               #########################################
########################################################################################################################

generate_initials_simple_thph <- function(def_sigma = 50, num = 10, nClass = 1, especificas = NULL,
                                     guess_PH = 3, guess_theta_w = 0, guess_certainty = 3, covariates = NULL) {
  if (nClass == 1) {
    testvals <- data.frame(
      PH = runif(num,0, 1),
      theta_w  = runif(num,-1, 1),
      sigma = rep(def_sigma, num))
    if (!is.null(covariates)) {
      testvals[, paste0("PH_",covariates)] <- 0
      testvals[, paste0("theta_w_",covariates)] <- 0 }
  } else {
    testvals <- data.frame( PH_1 = rep(0, num) )
    for (s in 1:nClass ) {
      testvals[,paste0("PH_",s)]     <- guess_PH[s] + runif(num, -guess_certainty, guess_certainty)
      testvals[testvals[,paste0("PH_",s)] < 0,paste0("PH_",s)]  <- 0.01
      testvals[,paste0("theta_w_",s)] <- guess_theta_w[s]  + runif(num, -guess_certainty, guess_certainty)
      testvals[,paste0("sigma_",s)]  <- rep(def_sigma, num)
      testvals[,paste0(especificas,"_", s)] <- matrix(runif(num*length(especificas), -1, 1), nrow = num, ncol = length(especificas))
      if (!is.null(covariates)) {
        testvals[, paste0("PH_",covariates,"_",s)] <- 0
        testvals[, paste0("theta_w_",covariates,"_",s)] <- 0 }
      }
    testvals[,paste0(especificas,"_", 1)] <- 0
  }
  return(testvals)
}

generate_initials_multi_thph <- function(def_sigma = 50, num = 10, nClass = 1, especificas = NULL,
                                     guess_PH = 0.5, guess_theta_w = 0, guess_certainty = 0.5, corners=F) {
  TH <- 1
  if (nClass == 1) {
    testvals <- data.frame(
      PH = runif(num,0, 1),
      theta_w  = runif(num,-1, 1))
    if (ntimes >= 2) {
      testvals[,mod_thetas] <- matrix(runif(num*length(mod_thetas), 0.1, 0.9), nrow = num, ncol = length(mod_thetas))
      sumthetas <- rowSums(testvals[,mod_thetas]) + runif(num,0.1, 0.9)
      testvals[,mod_thetas] <- sweep(testvals[,mod_thetas], 1, sumthetas, "/" ) }
    if (nexpenditures >= 1) {
      testvals[,mod_phis] <- matrix(runif(num*length(mod_phis), 0.1, 0.9), nrow = num, ncol = length(mod_phis))
      sumphis <- rowSums(testvals[,mod_phis]) + runif(num,0.1, 0.9)
      testvals[,mod_phis] <- sweep(testvals[,mod_phis], 1, sumphis, "/" )
      testvals[,mod_phis] <- sweep(testvals[,mod_phis], 1, testvals[,"PH"], "*" ) }
    testvals[, mod_sigmas] <- def_sigma
    testvals[,mod_rhos]   <- 0
  } else {
    testvals <- data.frame(
        PH_1 = rep(0, num))
    for (s in 1:nClass ) {
      if (corners) notdoing <- !act_mix[[s]]
      testvals[, paste0("PH_",s)]       <- guess_PH[s] + runif(num, -guess_certainty, guess_certainty)
      testvals[testvals[,paste0("PH_",s)] < 0,paste0("PH_",s)]  <- 0.01
      testvals[, paste0("theta_w_",s)]  <- guess_theta_w[s]  + runif(num, -guess_certainty, guess_certainty)
      if (ntimes >= 2) {
        testvals[,paste0(mod_thetas,"_",s)] <- matrix(runif(num*length(mod_thetas), 0.1, 0.9), nrow = num, ncol = length(mod_thetas))
        if (corners) if (sum(notdoing) > 0) {testvals[,paste0(mod_thetas,"_",s)][,notdoing[2:ntimes]] <- 0}
        sumthetas <- rowSums(testvals[,paste0(mod_thetas,"_",s)]) + runif(num,0.1, 0.9)
        testvals[,paste0(mod_thetas,"_",s)] <- sweep(testvals[,paste0(mod_thetas,"_",s)], 1, sumthetas, "/" ) }
      if (nexpenditures >= 1) {
        testvals[,paste0(mod_phis,"_",s)] <- matrix(runif(num*length(mod_phis), 0.1, 0.9), nrow = num, ncol = length(mod_phis))
        sumphis <- rowSums(testvals[,paste0(mod_phis,"_",s)]) + runif(num,0.1, 0.9)
        testvals[,paste0(mod_phis,"_",s)] <- sweep(testvals[,paste0(mod_phis,"_",s)], 1, sumphis, "/" )
        testvals[,paste0(mod_phis,"_",s)] <- sweep(testvals[,paste0(mod_phis,"_",s)], 1, testvals[,paste0("PH_",s)], "*" ) }
      testvals[, paste0(mod_sigmas, "_", s)] <- def_sigma
      testvals[,paste0(mod_rhos,"_",s)]   <- 0
      testvals[,paste0(especificas,"_", s)] <- matrix(runif(num*length(especificas), -1, 1), nrow = num, ncol = length(especificas)) }
    testvals[,paste0(especificas,"_", 1)] <- 0
  }
  return(testvals)
}

generate_initials_simple_tc_thph <- function(def_sigma = 50, num = 10, especificas = NULL, nClass=2,
                                     guess_PH = 3, guess_theta_w = 0, guess_certainty = 3, covariates = NULL) {

  testvals <- data.frame( PH_1 = rep(0, num) )
  for (s in 1:nClass ) {
    testvals[,paste0("PH_",s)]     <- guess_PH[s] + runif(num, -guess_certainty, guess_certainty)
    testvals[,paste0("theta_w_",s)] <- guess_theta_w[s]  + runif(num, -guess_certainty, guess_certainty)
    testvals[,paste0("sigma_",s)]  <- rep(def_sigma, num)
    testvals[,paste0(especificas,"_", s)] <- matrix(runif(num*length(especificas), -1, 1), nrow = num, ncol = length(especificas))
    if (!is.null(covariates)) {
      testvals[, paste0("PH_",covariates,"_",s)] <- 0
      testvals[, paste0("theta_w_",covariates,"_",s)] <- 0 }
    }
  testvals[,paste0(especificas,"_", 1)] <- 0
  return(testvals)
}

generate_initials_quadratic <- function(def_sigma = 10, num =10, nClass = 1, especificas = NULL, covariates = NULL, guess_certainty = 1) {
  vect <- c("alpha", "beta_w1", "beta_w2", "beta_ta1", "beta_ta2", "beta_Ec1", "beta_Ec2", "beta_wta", "beta_wEc", "beta_taEc")
  if (nClass == 1) {
    testvals <- as.data.frame(matrix( runif(num*length(vect), -guess_certainty, guess_certainty)   ,nrow=num, ncol= length(vect)))
    colnames(testvals) <- vect
    testvals[,"sigma"] <- def_sigma
    if (!is.null(covariates)) {
      testvals[, paste0("alpha_",covariates)] <- matrix(runif(num*length(covariates), -guess_certainty, guess_certainty), nrow = num, ncol=length(covariates))
      testvals[, paste0("beta_w1_",covariates)] <- matrix(runif(num*length(covariates), -guess_certainty, guess_certainty), nrow = num, ncol=length(covariates))
      testvals[, paste0("beta_ta1_",covariates)] <- matrix(runif(num*length(covariates), -guess_certainty, guess_certainty), nrow = num, ncol=length(covariates))
      testvals[, paste0("beta_Ec1_",covariates)] <- matrix(runif(num*length(covariates), -guess_certainty, guess_certainty), nrow = num, ncol=length(covariates))
    }
  } else {
    testvals <- as.data.frame(matrix( runif(num*length(vect)*nClass, -guess_certainty, guess_certainty) ,nrow=num, ncol= length(vect)*nClass))
    names(testvals) <- to_vec(for (s in 1:nClass) paste0(vect,"_",s) )
    testvals[,paste0("sigma_", 1:nClass)] <- def_sigma
    testvals[,to_vec(for (s in 1:nClass) paste0(especificas,"_", s))] <- matrix(runif(num*length(especificas)*nClass, -1, 1), nrow = num, ncol = length(especificas)*nClass)
    testvals[,paste0(especificas,"_", 1)] <- 0
  }
  return(testvals)
}

########################################################################################################################
###################################            REPORT VALUES OF TIME            ########################################
########################################################################################################################

report_values_of_time_thph <- function(best_model, dbs) {
  Tc  <- dbs[, "Tc"]
  Ec  <- dbs[, "Ec"]
  w   <- dbs[, "w"]
  ec  <- Ec / (w * (168 - Tc))

  TH             <- 1
  PH             <- as.numeric(best_model$estimate["PH"])
  theta_w        <- as.numeric(best_model$estimate["theta_w"])

  #### Valores del Tiempo
  thetaphiec <- PH + theta_w + (TH + theta_w)*ec
  auxsqrt <- sqrt(thetaphiec^2 - 4*theta_w*ec*(PH + TH + theta_w))
  Tw <- (168-Tc)*(thetaphiec + auxsqrt) / (2*(PH + TH + theta_w))

  dbs[, "cteVoLi"]  <- (w*Tw - Ec) / (168-Tw-Tc)
  dbs[, "cteVTAWi"] <- (w*Tw - Ec) / (Tw)
  cteVoL  <- mean(dbs$cteVoLi  , na.rm = T)
  cteVTAW <- mean(dbs$cteVTAWi , na.rm = T)
  coefVoL <-  TH / PH
  coefVTAW <- theta_w / PH

  dbs[, "VoL"]  <- coefVoL  * dbs[, "cteVoLi"]
  dbs[, "VTAW"] <- coefVTAW * dbs[, "cteVTAWi"]

  #### Intervalos de Confianza
  delta <- apollo_deltaMethod(best_model, list(expression=c(
    VoL      = glue('{cteVoL} *({TH} / PH)'),
    VTAW     = glue('{cteVTAW}*(theta_w / PH)'))))

  dbs[, "VoLse"]  <- VoLSE  <- delta[1,"s.e."]
  dbs[, "VTAWse"] <- VTAWSE <- delta[2,"s.e."]
  VoLmean  <- mean(dbs$VoL  )
  VTAWmean <- mean(dbs$VTAW )

  cat("VoL IC : [", round(VoLmean - 1.96*VoLSE, 4), ";", round(VoLmean + 1.96*VoLSE, 4),"]", end = "\n")
  cat("VTAW IC: [", round(VTAWmean - 1.96*VTAWSE, 4), ";", round(VTAWmean + 1.96*VTAWSE, 4),"]", end = "\n")
  var <- 1.96*sd(dbs$w)/sqrt(nrow(dbs))
  cat("observed wage rate:", mean(dbs$w), end = "\n")
  cat("observed wage rate IC: [",round(mean(dbs$w) - var, 4), ";", round(mean(dbs$w) + var, 4),"]", end = "\n")

  # dbs %>% ggplot(aes(x = VoL  / w, fill = "VoL")) + geom_histogram(bins = 100)
  # dbs %>% ggplot(aes(x = VTAW / w, fill = "VTAW")) + geom_histogram(bins = 100)
  # dbs <<- dbs
  return(delta)
}