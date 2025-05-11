pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", 
          "forcats", "lubridate", "haven", "reshape2","apollo")
invisible(lapply(pkgs, library, character.only=TRUE))

acts_corregidas <- c(
    # confiables:
    't_to',
    "t_tcnr_ce", # cuidados escenciales
    "t_tcnr_re", # cuidados relativos a la eseñanza
    "t_tcnr_oac", # otros cuidados
    #"t_tcnr_0_4",
    #"t_tcnr_5_14",
    #"t_tcnr_nna",  #(niños y adolescentes)
    #"t_tcnr_15_65",
    #"t_tcnr_66",
    #"t_tncr_psdf"
    "t_tdnr_psc",
    "t_tdnr_lv",
    "t_tdnr_lrc",
    "t_tdnr_mrm",
    "t_tdnr_admnhog",
    "t_tdnr_comphog",
    "t_tdnr_cmp",
    "t_tvaoh_tv", # voluntario a la comunidad
    "t_tvaoh_oh", #voluntario a otros hogares
    "t_cpaf_cp", # cuidados personales
    "t_cpag_comer",
    "t_cpag_dormir",
    "t_ed",
    "t_vsyo_csar", # convivencia social y actividades recreativas
    "t_vsyo_aa", # arte y aficiones

    't_mcm_leer',
    't_mcm_video',
    't_mcm_audio',
    't_mcm_computador',

    "t_tt1", # traslados enut 2015
    "t_tt2"  # traslados adicionales
  )

t_agregados <- c(
  "t_paid_work",
  "t_domestic_work",
  "t_care_work",
  "t_unpaid_voluntary",
  "t_education",
  "t_leisure",
  "t_personal_care",
  "t_meals",
  "t_sleep",
  "t_commute1",
  "t_commute2")

identificadores <- c("id_persona", "id_hog")
composicion_hogar <- c("parentesco",
                       "n_menores6",
                       "n_menores12",
                       "n_men15",
                       "n_menores18",
                       "n_menores",
                       "n_mayores",
                       "n_tiempo",
                       "n_trabajadores",
                       "n_profesionales",
                       "hay_tercera",
                       "n_personas",
                       "edad_promedio",
                       "tiene_hijos",
                       "en_pareja",
                       "vive_pareja")

sociodemograficas <- c("sexo",
                       "edad_anios",
                       "tramo_edad",
                       "NSE",
                       "nivel_escolaridad",
                       "estudia",
                       "trabaja",
                       "horas_trabajo",
                       "quintil",
                       "macrozona",
                       "region_ord",
                       "glosa_region",
                       "prop_hogar")
laborales <- c("cae", "cise", "ciuo_agrupada")
proveedores_externos <- c("servicio_domestico", "ayuda_cercanos", "fuentes_externas")
ingresos     <- c("ing_ocuppal", "ing_trab", "ing_jub_aps", "ing_g", "ing_t_hogar", "ing_t_pc", "ing_gpp", "ing_personal", "ingreso_hogar" ,"income_person_week")
tipo_muestra <- c("es_trabajador", "es_familia")

new_variables_prefilter<- function(data) { 
  # OJO, DEBEN ENTRAR TODOS INDEPENDIENTE DE SI DECLARAN TIEMPO O NO
  act_s_domestico <- c("a6","a14","a24","a32", "a1", "a19")
  recibe_servicio_domestico <- c("a1", "a19")
  
  # INFORMACIÓN DEL INDIVIDUO
  data <- data %>%
    mutate(
      menor_edad        = case_when(edad < 18 ~ 1, TRUE ~ 0),
      menor6            = case_when(edad < 6 ~ 1, TRUE ~ 0),
      menor12           = case_when(edad >= 6 & edad < 12 ~ 1, TRUE ~ 0),
      men15             = case_when(edad < 15 ~ 1, TRUE ~ 0),
      menor18           = case_when(edad >= 12 & edad < 18 ~ 1, TRUE ~ 0),
      menor25           = case_when(edad <= 25 ~ 1, TRUE ~ 0),
      mayor_edad        = case_when(edad >= 18 ~ 1, TRUE ~ 0),
      tercera_edad      = case_when(edad >= 60 ~ 1, TRUE ~ 0),
      nivel_escolaridad = nivel_educ, 
      # (1: 0,1; 2: 2,3,4; 3: 5,6,7,9; 4:8, 5:10,11,12)
  
      # tiene_hijos = case_when(t11_1_3 < 6 ~ 1, TRUE ~ 0),
      # trabaja     = case_when(k11_1_1 == 1 ~ 1, k13_1_1 == 1 ~ 1, TRUE ~ 0),
      trabaja     = case_when(o1 == 1 ~ 1, TRUE ~ 0),
      teletrabaja = case_when(to2_v_ds == 1 | to2_v_fds == 1 ~ 1, TRUE ~0),
      jornada_laboral =jor_to, 
      en_pareja   = case_when(c7 <= 3 ~ 1, TRUE ~0),
      vive_pareja = case_when(c8a == 1 ~ 1, TRUE ~ 0),
      dia_semana  = case_when(
        cd1 == "Lunes" ~ 1,
        cd1 == "Martes" ~ 2,
        cd1 == "Miercoles" ~ 3,
        cd1 == "Jueves" ~ 4,
        T ~ 5,
      ),
      dia_fin_semana = case_when(
        cd14 == "Sabado" ~ 6,
        T ~ 7)) %>%
    group_by(id_hog) %>% 
      mutate(tiene_hijos = case_when(
        sum(pco == 3) > 0 & pco <=2 ~ 1, T ~ 0)) %>% 
      ungroup()


  # COMPOSICIÓN DEL HOGAR
  data <- data %>%
    group_by(id_hog) %>%
    mutate(n_menores = sum(menor_edad,na.rm=TRUE),
           n_menores6 = sum(menor6,na.rm=TRUE),
           n_menores12 = sum(menor12,na.rm=TRUE),
           n_men15 = sum(men15,na.rm=TRUE),
           n_menores18 = sum(menor18,na.rm=TRUE),
           n_menores25 = sum(menor25, na.rm=TRUE),
           n_mayores = sum(mayor_edad,na.rm=TRUE),
           n_tercera = sum(tercera_edad,na.rm=TRUE),
           n_tiempo = sum(tiempo, na.rm = T),
           n_trabajadores = sum(trabaja, na.rm = T),
           n_profesionales = sum(nivel_escolaridad >=4 ),
           edad_promedio = mean(edad)) %>%
    mutate(hay_tercera = case_when(n_tercera > 1 ~ 1 & tercera_edad == 0,  T ~ 0),
           n_personas = max(hhsize, na.rm =T)) %>%
  ungroup()

  # AYUDAS QUE RECIBE EL HOGAR
  data <- data %>%
    mutate_at(act_s_domestico, ~ replace_na(.,0)) %>%
    mutate_at(act_s_domestico, ~ifelse(. == 96, 0, .)) %>%
    group_by(id_hog) %>%
    mutate(
      servicio_domestico = abs(max(a1) - 2),
      ayuda_cercanos = abs(max(a19) - 2)) %>%
    mutate(fuentes_externas = case_when(
      servicio_domestico  == 1 | ayuda_cercanos == 1 ~ 1, T ~ 0)) %>%
    ungroup()

  # PREPROCESAMIENTO DE INGRESOS
  data <- data %>%
    # filter(tiempo == 1) %>% # será relevante desaserme de las personas que no declaran tiempo antes?
    mutate_at(c("ing_ocuppal", "ing_trab", "ing_jub_aps", "ing_t_hogar", "ing_t_pc"),
              ~replace_na(.,0)) %>%
      mutate(ingresos_propios =  ing_trab + ing_jub_aps) %>%
      group_by(id_hog) %>%
      mutate(ing_g = (sum(ing_t_hogar)-sum(ing_trab)+sum(ing_jub_aps)) * (ing_t_hogar>0)) %>%
      mutate(prop_hogar = ingresos_propios / sum(ingresos_propios)) %>%
      mutate(ingreso_hogar = sum(ingresos_propios) + sum(ing_g)) %>%
      mutate(ing_gpp = sum(ing_g) * prop_hogar) %>%
    mutate(income_person_week = (sum(ingresos_propios) + ing_gpp) / n()) %>% ungroup() %>%
      mutate(ingreso_personal = ing_trab + ing_jub_aps + ing_gpp) %>%
      group_by(id_hog) %>%
    mutate(prop_hogar = ingreso_personal / sum(ingreso_personal)) %>%
    ungroup() %>%
    arrange(ingreso_hogar) %>%
      mutate(percentil = cumsum(fe_ch) / sum(fe_ch)) %>%
    mutate(quintil = case_when(percentil<= 0.20000001 ~1,
                               percentil<= 0.40000001 ~2,
                               percentil<= 0.60000001 ~3,
                               percentil<= 0.80000001 ~4,
                               T ~ 5 )) %>% dplyr::select(-c(percentil, fe_ch)) %>%
    mutate_at(c("prop_hogar", "ingreso_hogar", "ing_gpp", "income_person_week", "ingreso_personal"), ~replace_na(.,0))

  return(data)
}

na_completion <- function(data) {
  act_traslados_2015 <- c(
    "t_tto", # traslado al trabajo en la ocupación
    "t_ted", # traslado al establecimiento educacional
    "cp7_t", # ida y vuelta a establecimiento de salud
    "cp10_t"
    )
  act_traslados_2023 <- c(
    "td14_t", # ida y vuelta a hacer tramites
    "td17_t",
    "td20_t", # ida y vuelta a hacer compras
    "td23_t",
    "tc12_t", # llevar y buscar al establecimiento
    "tc15_t",
    "tc21_t",  # llevarle al establecimiento de salud
    "tc25_t",
    "tc31_t",  # llevarle al lugar de trabajo
    "tc34_t"
  )
  act_traslados_2015_ds <- paste0(act_traslados_2015, "_ds")
  act_traslados_2015_fds <- paste0(act_traslados_2015, "_fds")
  act_traslados_2023_ds <- paste0(act_traslados_2023, "_ds")
  act_traslados_2023_fds <- paste0(act_traslados_2023, "_fds")
  vars_todos <- c("bs1")

  acts <- c(
    # confiables:
    't_to',
    'cp1_t',

    # no confiables:
    "t_tcnr_ce", # cuidados escenciales
    "t_tcnr_re", # cuidados relativos a la eseñanza
    "t_tcnr_oac", # otros cuidados
    #"t_tcnr_0_4",
    #"t_tcnr_5_14",
    #"t_tcnr_nna",  #(niños y adolescentes)
    #"t_tcnr_15_65",
    #"t_tcnr_66",
    #"t_tncr_psdf"
    "t_tdnr_psc",
    "t_tdnr_lv",
    "t_tdnr_lrc",
    "t_tdnr_mrm",
    "t_tdnr_admnhog",
    "t_tdnr_comphog",
    "t_tdnr_cmp",
    "t_tvaoh_tv", # voluntario a la comunidad
    "t_tvaoh_oh", #voluntario a otros hogares
    "t_cpaf_cp", # cuidados personales
    "cp3_t", "cp4_t", "cp5_t", # crear t_cpag_comer
    "t_ed",
    "t_vsyo_csar", # convivencia social y actividades recreativas
    "t_vsyo_aa", # arte y aficiones
    "vs7_t",
    "vs8_t",
    "vs9_t",
    "vs10_t"
  )

  acts_ds <- paste0(c(acts, act_traslados_2015, act_traslados_2023), "_ds")
  acts_fds <- paste0(c(acts, act_traslados_2015, act_traslados_2023), "_fds")

  data <- data  %>%
    # TODO reemplazar los 96 con 0
    filter(tiempo == 1) %>%
    mutate(across(
        c(edad, sexo), function(x) ifelse(is.na(x), 85, x))) %>%
    filter(if_all(c('cp1_t_ds', 'cp1_t_fds', 't_to_ds', 't_to_fds'), ~ !.x  %in% 96)) %>%
    #filter(if_all(n_linea_p:te_ayuda_cercanos, ~ !.x  %in% 96)) %>%
    #filter(l110_1_1 != 1 & l122_1_1 != 1) %>% # en cama enfermos
    #filter(l111_1_1 != 1 & l123_1_1 != 1) %>%# pasaron cosas fuera de lo común
    #filter(m14_1_1 != 2 & m14_2_1 != 2) %>% # huelga
    #filter(k42_1_1 != 6) %>%# no buscó trabajo por estar embarazada
    #filter(!(l15_1_1 == 2 & l16_1_1 == 2 & l17_1_1 == 2 & (k42_1_1 == 2 | k42_1_1 == 3 | k42_1_1 == 5))) %>% # los de vacaciones
    #filter(!(l117_1_1 == 2 & l118_1_1 == 2 & l119_1_1 == 2 & (k42_1_1 == 2 | k42_1_1 == 3 | k42_1_1 == 5))) %>%
    #filter(k15_1_1 != 1) %>%
    mutate_at(c(acts_ds, acts_fds), ~ifelse(is.na(.), 0, .)) %>%
    mutate_at(c(acts_ds, acts_fds), ~ifelse(. > 25, 0, .)) %>%
    mutate_at(c(acts_ds, acts_fds), ~replace_na(.,0)) %>%
    mutate_at(c("o22a", "o22b", "o22c"), ~ifelse(. == 96,0, .)) %>%
    mutate_at(c("o22a", "o22b", "o22c"), ~replace_na(.,0)) %>%
    filter(if_all(all_of(vars_todos), ~ !.x  %in% 85)) %>%
    mutate(
      t_cpag_comer_ds = dplyr::select(., all_of(paste0(c("cp3", "cp4", "cp5"), "_t_ds"))) %>% rowSums(),
      t_cpag_comer_fds = dplyr::select(., all_of(paste0(c("cp3", "cp4", "cp5"), "_t_fds"))) %>% rowSums(),
      t_cpag_dormir_ds = cp1_t_ds,
      t_cpag_dormir_fds = cp1_t_fds,
    ) %>%
    mutate(
      t_cpaf_af_ds = t_cpaf_af_ds - t_cpag_dormir_ds- t_cpag_comer_ds,
      t_cpaf_af_fds = t_cpaf_af_fds - t_cpag_dormir_fds- t_cpag_comer_fds,
      t_cpaf_cp_ds = t_cpaf_cp_ds - cp7_t_ds- cp10_t_ds,
      t_cpaf_cp_fds = t_cpaf_cp_fds - cp7_t_fds- cp10_t_fds,
      t_tdnr_admnhog_ds = t_tdnr_admnhog_ds - td14_t_ds - td17_t_ds,
      t_tdnr_admnhog_fds = t_tdnr_admnhog_fds - td14_t_fds - td17_t_fds,
      t_tdnr_comphog_ds = t_tdnr_comphog_ds -  td20_t_ds - td23_t_ds,
      t_tdnr_comphog_fds = t_tdnr_comphog_fds - td20_t_fds - td23_t_fds,
      t_tcnr_re_ds = t_tcnr_re_ds - tc12_t_ds - tc15_t_ds,
      t_tcnr_re_fds = t_tcnr_re_fds - tc12_t_fds - tc15_t_fds,
      t_tcnr_oac_ds = t_tcnr_oac_ds - tc21_t_ds - tc25_t_ds - tc31_t_ds - tc34_t_ds,
      t_tcnr_oac_fds = t_tcnr_oac_fds - tc21_t_fds - tc25_t_fds - tc31_t_fds - tc34_t_fds,
      t_tt1_ds = dplyr::select(., all_of(act_traslados_2015_ds)) %>% rowSums(na.rm = T),
      t_tt1_fds = dplyr::select(., all_of(act_traslados_2015_fds)) %>% rowSums(na.rm = T),
      t_tt2_ds = dplyr::select(., all_of(act_traslados_2023_ds)) %>% rowSums(na.rm = T),
      t_tt2_fds = dplyr::select(., all_of(act_traslados_2023_fds)) %>% rowSums(na.rm = T),
      t_mcm_leer_ds = dplyr::select(., all_of(c('vs7_t_ds'))) %>% rowSums(na.rm = T),
      t_mcm_video_ds = dplyr::select(., all_of(c('vs8_t_ds'))) %>% rowSums(na.rm = T),
      t_mcm_audio_ds = dplyr::select(., all_of(c('vs9_t_ds'))) %>% rowSums(na.rm = T),
      t_mcm_computador_ds = dplyr::select(., all_of(c('vs10_t_ds'))) %>% rowSums(na.rm = T),

      t_mcm_leer_fds = dplyr::select(., all_of(c('vs7_t_fds'))) %>% rowSums(na.rm = T),
      t_mcm_video_fds = dplyr::select(., all_of(c('vs8_t_fds'))) %>% rowSums(na.rm = T),
      t_mcm_audio_fds = dplyr::select(., all_of(c('vs9_t_fds'))) %>% rowSums(na.rm = T),
      t_mcm_computador_fds = dplyr::select(., all_of(c('vs10_t_fds'))) %>% rowSums(na.rm = T)
    ) %>%
    mutate(
      t_total_ds = dplyr::select(., all_of(paste0(acts_corregidas, "_ds"))) %>% rowSums(na.rm = T),
      t_total_fds = dplyr::select(., all_of(paste0(acts_corregidas, "_fds"))) %>% rowSums(na.rm = T),
    ) %>%
    mutate(sexo = abs(sexo - 1))

  return(data)
}

outlier_detection_Vallejo <- function(data)  {

  data["dias_normal"] <- TRUE
  temp <- data %>%
    # group_by(años_escolaridad, k11_1_1, rango_edad,sexo) %>%
    group_by(quintil, trabaja, tramo_edad,sexo) %>%
    mutate(xi_ds = mean(t_total_ds) , sd_ds = sd(t_total_ds),
           xi_fds = mean(t_total_fds), sd_fds = sd(t_total_fds)) %>%
    mutate(gamma_ds = sd_ds/xi_ds, gamma_fds = sd_fds/xi_fds) %>%
    mutate(lower_limit_ds = xi_ds - sd_ds + gamma_ds/2,
           upper_limit_ds = xi_ds + sd_ds,
           lower_limit_fds = xi_fds - sd_fds + gamma_fds/2,
           upper_limit_fds = xi_fds + sd_fds) %>%
    ungroup() %>%
    mutate(dias_normal = case_when(
      dias_normal == FALSE ~ FALSE,
      t_total_ds  < lower_limit_ds  | t_total_ds  > upper_limit_ds  ~ FALSE,
      t_total_fds < lower_limit_fds | t_total_fds > upper_limit_fds ~ FALSE,
      TRUE ~ TRUE)) %>%
    filter(dias_normal == TRUE) %>%

  return(data)
}

new_variables_postfilter <- function(data) {
    ingresos     <- c("ing_ocuppal", "ing_trab", "ing_jub_aps", "ing_g", "ing_t_hogar", "ing_t_pc",
                      "ing_gpp", "ing_personal", "ingreso_hogar", "income_person_week")

    data <- data %>%
      mutate(horas_trabajo = case_when(o22c > 0 ~ o22c, T ~ 0),
             macrozona = case_when(region_ord <= 5  ~ "norte",
                                   region_ord == 7 ~ "metropolitana",
                                   region_ord <= 11  ~ "centro",
                                   TRUE ~ "sur"),
             edad_anios = edad,
             tramo_edad = case_when(
               tr_edad == 1 ~ "12-24",
               tr_edad == 2 ~ "25-44",
               tr_edad == 3 ~ "45-65",
               T            ~ "66+"), # 12-24, 25-45, 46-65, 66+
             estudia = abs(abs(e4-1)-1),
             n_menores  = case_when(n_menores  <= 3 ~ n_menores , T ~ 4),
             n_mayores = case_when(n_mayores <= 5 ~ n_mayores, T ~ 6),
             nivel_escolaridad = case_when( # todos los niveles son completos, excepto el inicial
               nivel_escolaridad == 1 ~ "ninguna",
               nivel_escolaridad == 2 ~ "primaria",
               nivel_escolaridad == 3 ~ "secundaria",
               nivel_escolaridad == 4 ~ "técnica",
               T                      ~ "universitaria"),
             parentesco = pco,
             ing_personal = ing_trab + ing_jub_aps + ing_gpp) %>%
      mutate_at(ingresos, ~. / 1000 / 4) %>% #/654
      mutate(ing_personal = round(ing_personal, 2))


    data <- data %>% mutate(
      cae = cae, # cae es igual a la anter
      cise = cise,
      ciuo_agrupada = ciuo_agrupada)

    cae_id <- c("Menor de 15 años"=1, "Ocupada(o)"=2, "Desocupada(o)"=3, "Inactiva(o)"=4, "Sin clasificar"=96)
    data$cae <- names(cae_id)[match(data$cae, cae_id)]

    return(data)
  }

data_to168hours <- function(data) {
  act_confiable     <- c("t_to", "t_cpag_dormir")
  act_no_confiable <- acts_corregidas[!acts_corregidas %in% act_confiable]

  data <- data %>%
    #mutate_at(c(act_confiable, act_no_confiable), ~ifelse(. <= 1/4,0 ,.)) %>%
    mutate(t_confiable =  dplyr::select(.,all_of(act_confiable)) %>% rowSums(na.rm = TRUE),
           t_no_confiable =  dplyr::select(.,all_of(act_no_confiable)) %>% rowSums(na.rm = TRUE)) %>%
    mutate(t_total = t_confiable + t_no_confiable)  %>%
    filter(t_confiable <= 168) %>%
    mutate(mult_no_confiable = (168- t_confiable)/(t_total - t_confiable)) %>%
    mutate_at(act_no_confiable, ~ .*mult_no_confiable) %>%
    mutate(
      t_total = dplyr::select(., all_of(c(act_confiable, act_no_confiable))) %>%
        rowSums(na.rm = TRUE))

  return(data)
}

impute_weekend <- function(data, twin_matrix) {

  acts <- list(
    "6" = paste0(acts_corregidas, "_sab"),
    "7" = paste0(acts_corregidas, "_dom"))

  finsemana <- paste0(acts_corregidas, "_fds")
  dias_semana <- paste0(acts_corregidas, "_ds")

  semana_completa <- acts_corregidas
  for (i in 6:7) {
    mask <- data$dia_fin_semana == i
    finde <- unlist(acts[as.character(i)])

    data[mask, finde]  <- data[mask, finsemana]
    suma <- rowSums(twin_matrix[!mask, mask])
    data[!mask, finde] <- as.matrix(twin_matrix[!mask, mask] ) %*% as.matrix(data[mask, finsemana])
    data[!mask, finde] <-  sweep(data[!mask, finde],1,suma,"/")
  }

  sabados <- unlist(acts["6"])
  domingos <- unlist(acts["7"])
  data_post <- data %>%
    #mutate_at(sabados , ~ifelse(. <= 1/12,0 ,.)) %>% # todo tiempo al que se le dedique menos de 5 minutos
    #mutate_at(domingos , ~ifelse(. <= 1/12,0 ,.)) %>%  # todo tiempo al que se le dedique menos de 5 minutos
    mutate(sum_sabados  = dplyr::select(., all_of(sabados)) %>% rowSums(na.rm = TRUE),
           sum_domingos = dplyr::select(., all_of(domingos)) %>% rowSums(na.rm = TRUE))%>%
    mutate_at(sabados , ~ . * 24 / sum_sabados)  %>%
    mutate_at(domingos, ~ . * 24 / sum_domingos) %>%
    mutate(sum_sabados = dplyr::select(., all_of(sabados)) %>% rowSums(na.rm = TRUE),
           sum_domingos = dplyr::select(., all_of(domingos)) %>% rowSums(na.rm = TRUE)) %>%
    mutate_at(dias_semana, ~case_when((o22b < 5) & (o22b > 0) ~ . *o22b, T ~ .* 5 ))

  temp <- data_post[,"t_to_ds"]
  data_post[,"t_to_ds"] <- 0
  data_post[,dias_semana] <- sweep(data_post[,dias_semana] , 1, rowSums(data_post[,dias_semana]), "/")
  data_post[,dias_semana] <- sweep(data_post[,dias_semana] , 1, as.numeric(unlist(24*5 - temp)), "*")
  data_post[,"t_to_ds"] <- temp

  data_post[,semana_completa] <- data_post[,dias_semana] + data_post[,sabados] + data_post[,domingos]
  data_post[,"t_total"] <-rowSums(data_post[,semana_completa])

  return(data_post)
}

diagnostico_trabajo <- function(data, plots, intercuartil) {
  diagnostico <- data %>%
    mutate(diferencia = t_to - o22c)

  cat("Número inicial de filas",nrow(diagnostico), end = "\n")

  if (plots) {
    ggplot(diagnostico, aes(x = o22a)) + geom_histogram(binwidth = 1)
    ggplot(diagnostico, aes(x = o22b)) + geom_histogram(binwidth = 1)
    ggplot(diagnostico, aes(x = o22c)) + geom_histogram(binwidth = 1)
    ggplot(diagnostico, aes(x = t_to)) + geom_histogram(binwidth = 1)
    ggplot(diagnostico, aes(x = diferencia)) + geom_histogram(binwidth = 1)
  }

  if (intercuartil) {
    cuantiles <- quantile(diagnostico[(diagnostico$t_to > 0) | (diagnostico$o22c > 0),]$diferencia, c(0.25, 0.75))
    intercuantil <- cuantiles[2] - cuantiles[1]
    diagnostico <- diagnostico %>% filter(diferencia < cuantiles[2] + intercuantil, diferencia > cuantiles[1] - intercuantil)
    cat("Número final de filas",nrow(diagnostico))
    return(diagnostico)
  } else {
    diagnostico <- diagnostico %>% filter(diferencia < 10, diferencia > -10)
    cat("Número final de filas",nrow(diagnostico))
    return(diagnostico)
  }
}

adjust_working_hours <- function(data, normalize = TRUE) {
  act_ajustar <- acts_corregidas[!acts_corregidas %in% c("t_to")]
  data <- data %>% mutate(t_to = horas_trabajo) # verificar con Astroza

  # if (normalize) {
  #   data <- data %>%
  #     mutate(t_ajustable = dplyr::select(.,act_ajustar) %>% rowSums(na.rm = TRUE)) %>%
  #     mutate(mult_ajuste= (168- tt)/(t_ajustable)) %>%
  #     mutate_at(act_ajustar, ~ .*mult_ajuste) %>%
  #     mutate(ttotal_ajuste = dplyr::select(., c("tt", act_ajustar))  %>% rowSums(na.rm = TRUE))
  # }

  return(data)
}

agregar_actividades <- function(data_post) {

  data_post <- data_post %>%
    mutate(
      t_paid_work = t_to, # free
      t_domestic_work = dplyr::select(., c("t_tdnr_psc", "t_tdnr_lv", "t_tdnr_lrc", "t_tdnr_mrm",
                                           "t_tdnr_admnhog", "t_tdnr_comphog", "t_tdnr_cmp",)) %>%
        rowSums(na.rm = TRUE), # commited/free
      t_care_work = dplyr::select(., c("t_tcnr_ce", "t_tcnr_re", "t_tcnr_oac",)) %>%
        rowSums(na.rm = TRUE),  # commited/free
      t_unpaid_voluntary = t_tvaoh_tv + t_tvaoh_oh, # free
      t_education = t_ed, # committed / free
      t_leisure = t_vsyo_csar + t_vsyo_aa + t_mcm_leer + t_mcm_video + t_mcm_audio + t_mcm_computador, # free
      t_personal_care = t_cpaf_cp, #committed/free
      t_meals         = t_cpag_comer, #committed/free
      t_sleep = t_cpag_dormir, #committed/free
      t_commute1 = t_tt1,
      t_commute2 = t_tt2)  # committed

  data_post <- data_post %>%
    mutate(
      es_trabajador = case_when(trabaja == 1 & cae == "Ocupada(o)" & t_to > 0 ~ 1, T ~ 0),
      es_familia = case_when(trabaja == 1 & cae == "Ocupada(o)" & t_to > 0 & vive_pareja == 1 & tiene_hijos ==1  ~ 1, T ~0))

  data25 <- data_post %>% dplyr::select(
    all_of(identificadores),
    all_of(tipo_muestra),
    dia_semana, dia_fin_semana,
    all_of(composicion_hogar),
    all_of(proveedores_externos),
    all_of(sociodemograficas),
    all_of(laborales),
    bs1:bs19,
    all_of(ingresos),
    all_of(acts_corregidas)) %>%
    mutate(t_total = dplyr::select(., all_of(acts_corregidas)) %>% rowSums(na.rm = TRUE))

  data25[,acts_corregidas]  = round(data25[,acts_corregidas],2)
  data25[,"temp"] = rowSums(data25[,acts_corregidas])
  data25[,"t_cpag_dormir"]       = data25[,"t_cpag_dormir"] - (data25[,"temp"] - 168)
  data25[,"temp"] = rowSums(data25[,acts_corregidas])
  data25 <- data25 %>% dplyr::select(-c("temp"))

  data11 <- data_post %>%
      dplyr::select(
        all_of(identificadores),
        all_of(tipo_muestra),
        dia_semana, dia_fin_semana,
        all_of(composicion_hogar),
        all_of(proveedores_externos),
        all_of(sociodemograficas),
        all_of(laborales),
        bs1:bs19,
        all_of(ingresos),
        all_of(t_agregados)) %>%
    mutate(t_total = dplyr::select(., all_of(t_agregados)) %>% rowSums(na.rm = TRUE))

  data11[,t_agregados]  = round(data11[,t_agregados],2)
  data11[,"temp"] = rowSums(data11[,t_agregados])
  data11[,"t_sleep"]  = data11[,"t_sleep"] - (data11[,"temp"] - 168)
  data11[,"temp"] = rowSums(data11[,t_agregados])

  return(list(data25 = data25, data11= data11))
}
###### MANTENER POR AHI ########
# grupos <- data %>%
#   group_by(dia_semana, años_escolaridad, trabaja, rango_edad,sexo) %>%
#   summarise(cuantos = n())
#######

imputacion_gastos <- function(data) {
  library("minpack.lm")
  data_hogar <- data %>%
    dplyr::select(id_hog, n_personas, n_menores6, n_menores12, n_menores18, n_menores, n_trabajadores, n_personas,
                  edad_promedio, quintil, macrozona, ingreso_hogar, income_person_week, n_profesionales) %>%
    distinct(id_hog, .keep_all = T) %>%
    mutate(n_personas_cut = case_when(n_personas >= 7 ~ 7, T ~ n_personas),
           n_menores6_cut = case_when(n_menores6 >= 3 ~ 3, T ~ n_menores12),
           n_menores12_cut = case_when(n_menores12 >= 2 ~ 3, T ~ n_menores12),
           n_menores18_cut = case_when(n_menores18 >= 2 ~ 3, T ~ n_menores18),
           n_trabajadores_cut = case_when(n_trabajadores >= 3 ~ 3, T ~ n_trabajadores),
           n_profesionales_cut = case_when(n_profesionales >= 2 ~ 3, T ~ n_profesionales))
  #hist(gastos$savings/gastos$ingresos, breaks =100)
  gastos  <<- read.csv("data/gastos.csv") %>% filter(ingreso_hogar > 0, savings/ingreso_hogar > -1)

  model   <<- apollo_loadModel("output/FMNL-epf-ix")
  lin_reg <<- lm(savings ~
                     ingreso_hogar +
                     (quintil==2) + (quintil==3) + (quintil==4) + (quintil==5) +
                     n_menores6_cut + n_menores12_cut + n_menores18_cut+ n_personas_cut +
                     n_trabajadores_cut  + n_profesionales_cut +
                     edad_promedio +
                     (macrozona == "norte") + (macrozona == "centro") + (macrozona == "sur"),
                 data = gastos)
  summary(lin_reg)
  codigos <<- c(g.alimentos = "01", g.vestimenta = "03",g.cuentas = "04", g.hogar = "05", g.salud = "06",
             g.transporte = "07", g.comunicaciones = "08", g.recreacion = "09", g.educacion = "10", g.restaurantes = "11")
  alts <<- c("alimentos", "vestimenta", "cuentas", "hogar", "salud", "transporte",
             "comunicaciones", "recreacion", "educacion", "restaurantes")

  mask <- data_hogar[, "ingreso_hogar"] >= 0
  data_hogar[mask,"savings"] <- predict(lin_reg, data_hogar[mask,])

  sum(mask)
  sum((data_hogar[mask,"savings"] < 0), na.rm = T )
  sum(( data_hogar[mask,"savings"] > data_hogar[mask,"ingreso_hogar"]), na.rm = T)
  data_hogar %>%
    distinct(id_hog, .keep_all = T) %>%
    group_by(quintil) %>%
    summarise(
      mean_savings = mean(savings, na.rm = TRUE),
      median_savings = median(savings, na.rm = TRUE),
      median_income = median(ingreso_hogar),
      min_hogar = min(ingreso_hogar))
  gastos %>%
    group_by(quintil) %>%
    summarise(
      mean_savings = mean(savings, na.rm = TRUE),
      median_savings = median(savings, na.rm = TRUE),
      median_income = median(ingreso_hogar),
      min_hogar = min(ingreso_hogar))

  data_hogar <-  data_hogar %>%
    mutate(savings = case_when(
      #savings < 0 ~ 0,
      savings > ingreso_hogar ~ ingreso_hogar, T ~ savings ))
  data_hogar$id_hogar <- 1:nrow(data_hogar)

  apollo_initialise()
  apollo_control <<- list(
    modelName       = "FMNL"   ,
    modelDescr      = "Fractional MNL model on time use data",
    indivID         = "id_hogar",
    outputDirectory = "output")
  apollo_probabilities <<- model$apollo_probabilities
  database <<- data_hogar
  database[,names(codigos)] <<- 0
  apollo_beta <<- model$apollo_beta
  apollo_fixed <<- model$apollo_fixed
  apollo_inputs <<- apollo_validateInputs()
  predictions_base <<- apollo_prediction(model, apollo_probabilities, apollo_inputs)
  data_hogar[,alts] = predictions_base[,alts]
  data_hogar[,alts] = sweep(data_hogar[,alts], 1, unlist(data_hogar[,"ingreso_hogar"] - data_hogar[,"savings"]), "*")

  data <- merge(data, data_hogar[,c("id_hog", "savings", alts)], by = "id_hog", all.x = TRUE)
  data[,c("savings", alts)] = sweep(data[,c("savings", alts)], 1, unlist(data[,"prop_hogar"]), "*")
  data[,"total_expenses"] = data[,"ing_personal"] - data["savings"]

  data[,alts] = round(data[,alts], 2)

  #ingresos     <- c("ing_ocuppal", "ing_trab", "ing_jub_aps", "ing_g", "ing_mon", "ing_mon_pc", "ing_gpp", "ing_personal", "income_person_week")
  # data[, ingresos] = data[, ingresos] * 1000
  # data[, alts] = data[, alts] * 1000
  # data[, "savings"] = data[, "savings"] * 1000
  # data[, "w"] = data[, "w"] * 1000

  return(data)
}
