rm(list = ls())  # .rs.restartR()
pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats", "lubridate", "haven", "reshape2","apollo", "comprehenr", "Hmisc")
est_set <- list(writeIter = FALSE, silent = T, maxIterations=500, scaleHessian = F, scaleAfterConvergence = F, estimationRoutine = "bfgs",
                bgw_settings = list(maxFunctionEvals = 1000), validateGrad  = FALSE)
invisible(lapply(pkgs, library, character.only=TRUE))
library(lavaan)
library(glue)
set.seed(42)

model_data <- haven::read_dta("../data/enut-ii-11G.dta") %>%
  filter(edad_años >= 18, es_trabajador == 1, ing_personal > 0, w > 0, total_expenses/ing_personal <= 5) %>%
  mutate(ta = 168,
         Tw = t_paid_work,
         Tc = 168 - t_leisure -  t_paid_work,
         Ec = cuentas + hogar + salud + transporte + educacion - ing_gpp - ing_jub_aps,
         recreation_expenses = recreacion , # + restaurantes + comunicaciones + vestimenta + alimentos
         t_commute = t_commute1 + t_commute2)  %>%
  mutate(ec = Ec / (w * (ta-Tc)))
model_data = model_data %>% mutate(wbar = mean(w), ta= 168 - Tc, tabar = mean(168-Tc), Ecbar = mean(Ec))
model_data <- type.convert(model_data, as.is =TRUE)

model_data <- model_data %>%
  mutate(edad1 = (tramo_edad == "12-24"),
         edad2 = (tramo_edad == "25-44"),
         edad3 = (tramo_edad == "45-65"),
         edad4 = (tramo_edad == "66+"),
         q1    = (quintil == 1),
         q2    = (quintil == 2),
         q3    = (quintil == 3),
         q4    = (quintil == 4),
         q5    = (quintil == 5),
         norte = (macrozona == "norte"),
         centro = (macrozona == "centro"),
         sur = (macrozona == "sur"),
         female = (sexo==1),
         edad = edad_años
  )

sociodemograficas <- 'edad1 + edad3 + edad4 + q2 + q3 + q4 + q5'
regionales <- 'norte + centro + sur'
formulation <- glue('t_sleep    ~  {sociodemograficas}+ {regionales}
                  t_commute       ~  t_sleep + {sociodemograficas}+ {regionales}
                  t_paid_work     ~  t_sleep + t_commute + {sociodemograficas}+ {regionales}
                  t_education     ~  t_sleep + t_commute + t_paid_work + {sociodemograficas} + {regionales}
                  t_meals         ~  t_sleep + t_commute + t_paid_work + t_education + {sociodemograficas} + {regionales}
                  t_care_work     ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + {sociodemograficas} + {regionales}
                  t_domestic_work ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + t_care_work + {sociodemograficas} + {regionales}
                  t_leisure       ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + t_care_work + t_domestic_work + {sociodemograficas} + {regionales} ')

fit = sem(formulation, data = model_data, estimator = "MLM",   meanstructure = T)
summary(fit)

sociodemograficas <- 'female + edad1 + edad3 + edad4 + q2 + q3 + q4 + q5'
regionales <- 'norte + centro + sur'
formulation <- glue('t_sleep    ~  {sociodemograficas}+ {regionales}
                  t_commute       ~  t_sleep + {sociodemograficas}+ {regionales}
                  t_paid_work     ~  t_sleep + t_commute + {sociodemograficas}+ {regionales}
                  t_education     ~  t_sleep + t_commute + t_paid_work + {sociodemograficas} + {regionales}
                  t_meals         ~  t_sleep + t_commute + t_paid_work + t_education + {sociodemograficas} + {regionales}
                  t_care_work     ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + {sociodemograficas} + {regionales}
                  t_domestic_work ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + t_care_work + {sociodemograficas} + {regionales}
                  t_leisure       ~  t_sleep + t_commute + t_paid_work + t_education + t_meals + t_care_work + t_domestic_work + {sociodemograficas} + {regionales} ')

fit = sem(formulation, data = model_data, estimator = "MLM",   meanstructure = T)
summary(fit)

model_data <- haven::read_dta("../data/enut-ii-22G.dta") %>%
  filter(edad_años >= 18, es_trabajador == 1, ing_personal > 0, w > 0, total_expenses/ing_personal <= 5) %>%
  mutate(ta = 168,
         Tw = t_to,
         Tl = t_vsyo_aa + t_vsyo_csar,
         Tc = 168 - t_vsyo_aa + t_vsyo_csar - t_mcm - t_to,
         Ec = cuentas + hogar + salud + transporte + educacion - ing_gpp - ing_jub_aps,
         recreation_expenses = recreacion +  restaurantes) %>% #+ comunicaciones + vestimenta + alimentos)  %>%
  mutate(ec = Ec / (w * (ta-Tc)))
model_data = model_data %>% mutate(wbar = mean(w), ta= 168 - Tc, tabar = mean(168-Tc), Ecbar = mean(Ec))
model_data <- type.convert(model_data, as.is =TRUE)

model_data <- model_data %>%
  mutate(edad1 = (tramo_edad == "12-24"),
         edad2 = (tramo_edad == "25-44"),
         edad3 = (tramo_edad == "45-65"),
         edad4 = (tramo_edad == "66+"),
         q1    = (quintil == 1),
         q2    = (quintil == 2),
         q3    = (quintil == 3),
         q4    = (quintil == 4),
         q5    = (quintil == 5),
         norte = (macrozona == "norte"),
         centro = (macrozona == "centro"),
         sur = (macrozona == "sur"),
         female = (sexo==1),
         edad = edad_años
  ) %>%
  mutate(
    Tl_norte = Tl * norte,
    Tl_centro = Tl * centro,
    Tl_sur = Tl * sur,

  )

formulation <- glue('
  Tw                  ~ {sociodemograficas}+ Tc
  Tl                  ~ Tw + {sociodemograficas} + Tc
  recreation_expenses ~ Tw + Tl + Tl_norte + Tl_centro + Tl_sur + {sociodemograficas} + Tc + Ec + ing_personal
  ')
# OBS: t_leisure  no puede relacionarse con Tc, dado que este ya se encuentra incluido en t_paid_work
fit = sem(formulation, data = model_data, estimator = "ML",   meanstructure = T)
summary(fit)