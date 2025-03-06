rm(list = ls())  # .rs.restartR()
pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", 
          "stringr", "forcats", "lubridate", "haven", "reshape2")
invisible(lapply(pkgs, library, character.only=TRUE))
options(dplyr.summarise.inform = FALSE)

source("data_processing/processing_functions.R")

## PREPROCESAMIENTO DE DATOS
data <- haven::read_dta("data/raw/ii-enut-bdd-stata.dta")
data <- new_variables_prefilter(data)
data <- na_completion(data)
data <- new_variables_postfilter(data)

# TODO darle relevancia a nivel socioeconomico de la unidad de muestreo NSE

grupos <- data %>%
  group_by(quintil, trabaja, tramo_edad, sexo) %>%
  summarise(cuantos = n())

# library(fitdistrplus)
# descdist(data$ttotal_ds_original, boot = 1000)
# fit_gamma  <- fitdist(data$ttotal_ds_original, "gamma")
# fit_normal <- fitdist(data$ttotal_ds_original, "norm")
# fit_lognormal <- fitdist(data$ttotal_ds_original, "lnorm")
#
# par(mfrow = c(2, 2))
# plot.legend <- c("Gamma","Normal","Lognormal")
# denscomp(list(fit_gamma,fit_normal,fit_lognormal), legendtext = plot.legend)
# qqcomp(list(fit_gamma,fit_normal,fit_lognormal), legendtext = plot.legend)
# cdfcomp(list(fit_gamma,fit_normal,fit_lognormal), legendtext = plot.legend)
# ppcomp(list(fit_gamma,fit_normal,fit_lognormal), legendtext = plot.legend)
# par(mfrow=c(1,1))

ggplot(data, aes(x = t_total_ds)) + geom_histogram(binwidth = 1)
ggplot(data, aes(x = t_total_fds)) + geom_histogram(binwidth = 1)

data_outliers <- outlier_detection_Vallejo(data)
#data_outliers <- outlier_detection_mahalanobis(data)
#data_outliers <- outlier_detection_quantiles(data)

ggplot(data_outliers, aes(x = t_total_ds)) + geom_histogram(binwidth = 1)
ggplot(data_outliers, aes(x = t_total_ds)) + geom_histogram(binwidth = 1)

data_outliers <- data_outliers %>% dplyr::select(
  all_of(identificadores),
  dia_semana,
  dia_fin_semana,
  o22a, o22b, o22c,
  #c5:c14,
  all_of(composicion_hogar),
  all_of(proveedores_externos),
  all_of(sociodemograficas),
  all_of(laborales),
  all_of(ingresos),
  all_of(paste0(acts_corregidas, "_ds")),
  all_of(paste0(acts_corregidas, "_fds")),
  bs1:bs19)

# ------------------------------------------------------------------------------
haven::write_dta(data_outliers %>% arrange(id_persona), "data/raw/ENUT_PRE_WEEKEND_IMPUTATION.dta")
write_csv(data_outliers %>% arrange(id_persona), "data/raw/ENUT_PRE_WEEKEND_IMPUTATION.csv")
data <- haven::read_dta("data/raw/ENUT_PRE_WEEKEND_IMPUTATION.dta") %>% arrange(id_persona)
# ------------------------------------------------------------------------------

twin_matrix <- read_csv("data/raw/matriz_gemelos.csv.gzip", col_names = F)
# twin_matrix <- t(as.matrix(twin_matrix))
data_post <- impute_weekend(data, twin_matrix) # TODO arreglar esta función para la nueva implementación de datos
rm(twin_matrix)
data_post <- diagnostico_trabajo(data_post, F, T)
data_post <- adjust_working_hours(data_post)
data_post <- data_to168hours(data_post) # TODO probar distintos posicionamientos

### Continuar con la imputación de gastos, está la opción de usar fractional logit para imputar los gastos
#aux <- data_post %>% filter(ing_trab > 0, t_to > 0, trabaja == 1) %>%
#  mutate(w = ing_trab / t_to) %>%
#  dplyr::select(ing_personal, ing_trab, w, tt, edad_años, trabaja)

data_descargable <- agregar_actividades(data_post)
data22 <- data_descargable[["data22"]] %>% mutate(w = ing_trab / t_to)
data11 <- data_descargable[["data11"]] %>% mutate(w = ing_trab / t_paid_work)
haven::write_dta(data22, "data/enut-ii-22.dta")
haven::write_dta(data11, "data/enut-ii-11.dta")
write_csv(data22, "data/enut-ii-22.csv")
write_csv(data11, "data/enut-ii-11.csv")

source("data_processing/processing_functions.R")

data22 = haven::read_dta( "data/enut-ii-22.dta")
data11 = haven::read_dta( "data/enut-ii-11.dta")

table(data22$es_trabajador)
table(data22$es_familia)

## ----- Análisis datos resultantes -----------------
data11G <- imputacion_gastos(data11)
data22G <- imputacion_gastos(data22)

haven::write_dta(data22G, "data/enut-ii-22G.dta")
haven::write_dta(data11G, "data/enut-ii-11G.dta")
write_csv(data22G, "data/enut-ii-22G.csv")
write_csv(data11G, "data/enut-ii-11G.csv")

ggplot(data11G, aes(x = w)) + geom_histogram(bins = 50)

#library(fitdistrplus)
#fg <- fitdist(data$t_total_ds, "gamma")
#ks.test(data$t_total_ds, "pgamma", shape = fg$estimate["shape"], rate = fg$estimate["rate"])
# No sigue gamma

#fn <- fitdist(data$t_total_ds, "norm")
#ks.test(data$t_total_ds, "pnorm", mean = fg$estimate["mean"], sd = fg$estimate["sd"])

#fl <- fitdist(data$t_total_ds, "lnorm")
#ks.test(data$t_total_ds, "plnorm", meanlog = fg$estimate["meanlog"], sdlog = fg$estimate["sdlog"])

#fw <- fitdist(data$t_total_ds, "weibull")
#ks.test(data$t_total_ds, "pweibull", shape = fg$estimate["shape"], scale = fg$estimate["scale"])

#plot.legend <- c("gamma", "norm", "lnorm", "weibull")
#denscomp(list(fg, fn, fl, fw), legendtext = plot.legend)
#qqcomp(list(fg, fn, fl, fw), legendtext = plot.legend)
#cdfcomp(list(fg, fn, fl, fw), legendtext = plot.legend)
#ppcomp(list(fg, fn, fl, fw), legendtext = plot.legend)
# https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
#descdist(data$t_total_ds, discrete = FALSE)

# voy aqui xdxd -----------
