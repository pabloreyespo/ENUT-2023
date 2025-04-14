rm(list = ls())  # .rs.restartR()
pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats", "lubridate", "haven", "reshape2","apollo")
invisible(lapply(pkgs, library, character.only=TRUE))

codigos <- c(g.alimentos = "01", g.vestimenta = "03",g.cuentas = "04", g.hogar = "05", g.salud = "06",
             g.transporte = "07", g.comunicaciones = "08", g.recreacion = "09", g.educacion = "10", g.restaurantes = "11")

vector_gastos <- c("g.alimentos","g.vestimenta","g.cuentas","g.hogar",
                   "g.salud","g.transporte","g.comunicaciones","g.recreacion","g.educacion",
                   "g.restaurantes","total_expenses")

data_ingresos <- haven::read_dta("data/raw/base-personas-ix-epf-stata.dta") 
data_gastos <- haven::read_dta("data/raw/base-gastos-ix-epf-stata.dta")

######################################################################################################
######################                      INGRESOS                          ######################## 
######################################################################################################

# codigos especialeS: -66, -77, (-99, missing)
data_ingresos <- data_ingresos %>% # sprincipal: jefe de hogar
  dplyr::select(folio, 
                fe, 
                macrozona, 
                npersonas, 
                sprincipal, 
                parentesco, 
                sexo, 
                edad, 
                ecivil, 
                ecompras, 
                edue, 
                eduactual, 
                educurso, 
                edunivel, 
                edutermina, 
                edudependencia, 
                cae, 
                caeg, 
                ing_total_hogar_hd, 
                salud, 
                contribuciones, 
                prevision, 
                impuesto, 
                ot05, 
                ot06, 
                gastot_hd, 
                gastot_hd_ai, 
                ing_disp_td, 
                ing_disp_ti, 
                #ing_disp_asal_hd, 
                #ing_disp_hon_hd, 
                #ing_neto_empleador_hd, 
                #ing_disp_cpneg_hd,
                # falta INGDNP_HD
                #ing_disp_prof_hd,
                #ing_disp_asal_hd, #INGDOTA
                #ing_disp_jub_hd,
                #ing_disp_ot_ind,
                #ing_propiedad, 
                #ing_financieros,
                #transferencias_emitidas,
                #transferencias_recibidas,
                ing_disp_hog_hd,
                ing_disp_hd) %>%
  #mutate(ing_disp_prof_hd = case_when(ing_disp_prof_hd < 0 ~ 0, T ~ ing_disp_prof_hd)) %>%
  mutate(ot05 = case_when(ot05 < 0 ~ 0, T ~ ot05), ot06 = case_when(ot06 < 0 ~ 0, T ~ ot06),
         edad = case_when(edad < 0 ~ 0, T ~ edad)) %>%
  mutate(parentesco = case_when(parentesco == 1 ~ 1, T ~ 0)) %>% 
  filter_all(all_vars(. != -99)) %>%
  filter_all(all_vars(. != -88)) %>%
  arrange(ing_disp_hog_hd) %>%
  mutate(percentil = cumsum(fe) / sum(fe)) %>%
  mutate(quintil = case_when(percentil<= 0.20000001 ~1,
                             percentil<= 0.40000001 ~2,
                             percentil<= 0.60000001 ~3,
                             percentil<= 0.80000001 ~4,
                             T ~ 5 )) %>% dplyr::select(-c(percentil, fe)) %>%
  group_by(folio) %>%
  summarise(n_menores6             = sum( edad < 6 ),
            n_personas             = max(npersonas),
            macrozona              = max(macrozona),
            n_menores12            = sum( edad >= 6  & edad < 12),
            n_menores18            = sum( edad >= 12 & edad < 18),
            n_trabajadores         = sum( (ing_disp_td + ing_disp_ti) > 0),
            income_person_week     = mean(ing_disp_hog_hd) / n(),
            edad_promedio          = mean(edad),
            n_profesionales        = sum(edunivel >= 14 & !(edunivel == 14 & edutermina == 2)),
            ingreso_hogar          = sum(ing_disp_hd),
            quintil                = max(quintil),
            ing_disp_hog_hd        = max(ing_disp_hog_hd)) %>%
  mutate(macrozona = case_when(
    macrozona == 1 ~ "norte",
    macrozona == 2 ~ "metropolitana",
    macrozona == 3 ~ "centro",
    macrozona == 4 ~ "sur")) %>%
  ungroup() %>% distinct(folio,.keep_all = TRUE)


data_gastos <- data_gastos %>%
  filter(ccif!="04.2.1.01.01" & ccif!="04.2.2.01.01" & ccif!="04.2.2.01.02" & folio!="") %>%
  drop_na(folio) %>%
  mutate(codigo = d ) %>% #paste(D,G, sep =".")) %>%
  mutate(codigo = case_when(
    d=="09" & g=="1" ~ "08",
    d=="02"          ~ "01", # alcoho
    d=="12" & g=="1" ~ "06",
    d=="12" & g=="2" ~ "09",
    d=="12" & g=="3" ~ "03",
    d=="12" & g=="4" ~ "04",
    d=="12" & g=="5" ~ "04",
    d=="12" & g=="6" ~ "04",
    d=="12" & g=="7" ~ "04",
    d=="12" & g=="8" ~ "04",
    T ~ codigo )) %>%
  group_by(folio, codigo) %>%
  summarise(monto_total = sum(gasto)) %>%
  mutate(codigo = names(codigos)[match(codigo, codigos)]) %>%
  ungroup() %>%
  dcast(folio ~ codigo, value.var="monto_total", fun.aggregate = sum) %>% dplyr::select(-c("NA")) %>%
  mutate_at(names(codigos), ~replace_na(.,0))

data_ambos <- merge(data_ingresos, data_gastos, "folio") %>%
  arrange(folio) %>%
  mutate(total_expenses = dplyr::select(.,names(codigos)) %>% rowSums()) %>% # legal.discounts
  mutate(savings = ingreso_hogar - total_expenses) %>%
filter(total_expenses > 0 & ingreso_hogar > 0 ) # & EDAD > 0 & EDUE > 0 & EDUNIVEL>0

#data_ambos <- data_ambos %>% filter(JHOGAR == 1)
#data_ambos %>% group_by(quintil, macrozona) %>% summarise(minim = min(ing_disp_hog_hd)) %>% arrange(quintil, macrozona)
#data_ambos %>% group_by(quintil, macrozona) %>% summarise(minim = min(ingreso_hogar)) %>% arrange(quintil, macrozona)

gastos <- data_ambos %>% mutate(id_hogar = 1:n()) %>%
  dplyr::select(folio, id_hogar, all_of(vector_gastos), savings, n_personas, 
                macrozona, quintil,n_menores6,
                n_menores12, n_menores18, n_trabajadores,edad_promedio, 
                n_profesionales, income_person_week, ingreso_hogar)

gastos[,c(vector_gastos,"ingreso_hogar", "savings","income_person_week" )] <- 
  gastos[,c(vector_gastos,"ingreso_hogar", "savings","income_person_week" )] / 4000 * 1.18 # inflación entre periodos de muestreo1.025 # (4 semanas, miles de pesos, 1.025 la inflación)

gastos <- gastos %>%
  mutate(n_personas_cut = case_when(n_personas >= 7 ~ 7, T ~ n_personas),
         n_menores6_cut = case_when(n_menores6 >= 3 ~ 3, T ~ n_menores12),
         n_menores12_cut = case_when(n_menores12 >= 2 ~ 3, T ~ n_menores12),
         n_menores18_cut = case_when(n_menores18 >= 2 ~ 3, T ~ n_menores18),
         n_trabajadores_cut = case_when(n_trabajadores >= 3 ~ 3, T ~ n_trabajadores),
         n_profesionales_cut = case_when(n_profesionales >= 2 ~ 3, T ~ n_profesionales))
write.csv(gastos, "data/gastos.csv")

########################################################################################################################
################                                 MODELOS                                                ################
########################################################################################################################

alts <- c("alimentos", "vestimenta", "cuentas", "hogar", "salud", "transporte",
          "comunicaciones", "recreacion", "educacion", "restaurantes")

temp = sweep(gastos[,paste0("g.", alts)],1,rowSums(gastos[,paste0("g.", alts)]),"/")
colMeans(temp)
apollo_initialise()
apollo_control <- list(
  modelName       = "FMNL-epf-ix"   ,
  modelDescr      = "Fractional MNL model on time use data",
  indivID         = "id_hogar", 
  outputDirectory = "output", 
  nCores = 44)

database <- gastos
database[,names(codigos)] <- database[, names(codigos)] / rowSums(database[, names(codigos)])

apollo_beta <- c(asc_alimentos          = 0, bnpersonas_alimentos          = 0, bnmenores6_alimentos = 0          , bnmenores12_alimentos = 0          , bntrabajadores_alimentos = 0         , bnprofesionales_alimentos = 0          , bedadpromedio_alimentos = 0          , bq2_alimentos = 0          , bq3_alimentos = 0          , bq4_alimentos = 0          , bq5_alimentos = 0          , bnorte_alimentos = 0          , bcentro_alimentos = 0          , bsur_alimentos = 0          ,
                 asc_vestimenta         = 0, bnpersonas_vestimenta         = 0, bnmenores6_vestimenta = 0         , bnmenores12_vestimenta = 0         , bntrabajadores_vestimenta = 0        , bnprofesionales_vestimenta = 0         , bedadpromedio_vestimenta = 0         , bq2_vestimenta = 0         , bq3_vestimenta = 0         , bq4_vestimenta = 0         , bq5_vestimenta = 0         , bnorte_vestimenta = 0         , bcentro_vestimenta = 0         , bsur_vestimenta = 0         ,
                 asc_cuentas            = 0, bnpersonas_cuentas            = 0, bnmenores6_cuentas = 0            , bnmenores12_cuentas = 0            , bntrabajadores_cuentas = 0           , bnprofesionales_cuentas = 0            , bedadpromedio_cuentas = 0            , bq2_cuentas = 0            , bq3_cuentas = 0            , bq4_cuentas = 0            , bq5_cuentas = 0            , bnorte_cuentas = 0            , bcentro_cuentas = 0            , bsur_cuentas = 0            ,
                 asc_hogar              = 0, bnpersonas_hogar              = 0, bnmenores6_hogar = 0              , bnmenores12_hogar = 0              , bntrabajadores_hogar = 0             , bnprofesionales_hogar = 0              , bedadpromedio_hogar = 0              , bq2_hogar = 0              , bq3_hogar = 0              , bq4_hogar = 0              , bq5_hogar = 0              , bnorte_hogar = 0              , bcentro_hogar = 0              , bsur_hogar = 0              ,
                 asc_salud              = 0, bnpersonas_salud              = 0, bnmenores6_salud = 0              , bnmenores12_salud = 0              , bntrabajadores_salud = 0             , bnprofesionales_salud = 0              , bedadpromedio_salud = 0              , bq2_salud = 0              , bq3_salud = 0              , bq4_salud = 0              , bq5_salud = 0              , bnorte_salud = 0              , bcentro_salud = 0              , bsur_salud = 0              ,
                 asc_transporte         = 0, bnpersonas_transporte         = 0, bnmenores6_transporte = 0         , bnmenores12_transporte = 0         , bntrabajadores_transporte = 0        , bnprofesionales_transporte = 0         , bedadpromedio_transporte = 0         , bq2_transporte = 0         , bq3_transporte = 0         , bq4_transporte = 0         , bq5_transporte = 0         , bnorte_transporte = 0         , bcentro_transporte = 0         , bsur_transporte = 0         ,
                 asc_comunicaciones     = 0, bnpersonas_comunicaciones     = 0, bnmenores6_comunicaciones = 0     , bnmenores12_comunicaciones = 0     , bntrabajadores_comunicaciones = 0    , bnprofesionales_comunicaciones = 0     , bedadpromedio_comunicaciones = 0     , bq2_comunicaciones = 0     , bq3_comunicaciones = 0     , bq4_comunicaciones = 0     , bq5_comunicaciones = 0     , bnorte_comunicaciones = 0     , bcentro_comunicaciones = 0     , bsur_comunicaciones = 0     ,
                 asc_recreacion         = 0, bnpersonas_recreacion         = 0, bnmenores6_recreacion = 0         , bnmenores12_recreacion = 0         , bntrabajadores_recreacion = 0        , bnprofesionales_recreacion = 0         , bedadpromedio_recreacion = 0         , bq2_recreacion = 0         , bq3_recreacion = 0         , bq4_recreacion = 0         , bq5_recreacion = 0         , bnorte_recreacion = 0         , bcentro_recreacion = 0         , bsur_recreacion = 0         ,
                 asc_educacion          = 0, bnpersonas_educacion          = 0, bnmenores6_educacion = 0          , bnmenores12_educacion = 0          , bntrabajadores_educacion = 0         , bnprofesionales_educacion = 0          , bedadpromedio_educacion = 0          , bq2_educacion = 0          , bq3_educacion = 0          , bq4_educacion = 0          , bq5_educacion = 0          , bnorte_educacion = 0          , bcentro_educacion = 0          , bsur_educacion = 0          ,
                 asc_restaurantes       = 0, bnpersonas_restaurantes       = 0, bnmenores6_restaurantes = 0       , bnmenores12_restaurantes = 0       , bntrabajadores_restaurantes = 0      , bnprofesionales_restaurantes = 0       , bedadpromedio_restaurantes = 0       , bq2_restaurantes = 0       , bq3_restaurantes = 0       , bq4_restaurantes = 0       , bq5_restaurantes = 0       , bnorte_restaurantes = 0       , bcentro_restaurantes = 0       , bsur_restaurantes = 0       )

apollo_fixed <- c("asc_vestimenta", "bnpersonas_vestimenta", "bnmenores6_vestimenta", "bnmenores12_vestimenta", "bntrabajadores_vestimenta",
                  "bnprofesionales_vestimenta", "bedadpromedio_vestimenta", "bq2_vestimenta", "bq3_vestimenta", "bq4_vestimenta", "bq5_vestimenta", 
                  "bnorte_vestimenta", "bcentro_vestimenta", "bsur_vestimenta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  V <- list()
  for (alt in c("alimentos", "vestimenta", "cuentas", "hogar", "salud", "transporte","comunicaciones", "recreacion", "educacion", "restaurantes")) {
    V[[alt]] <- get(paste0("asc_", alt)) +
      get(paste0("bnpersonas_", alt))      * n_personas_cut +#(n_personas >= 4 )                   +
      get(paste0("bnmenores6_", alt))      * n_menores6_cut +#(n_menores12 >= 1)                  +
      get(paste0("bnmenores12_", alt))     * n_menores12_cut +#(n_menores12 >= 1)                  +
      get(paste0("bntrabajadores_", alt))  * n_trabajadores_cut +#(n_trabajadores >= 2 )               +
      get(paste0("bnprofesionales_", alt)) * n_profesionales_cut +#(n_trabajadores >= 2 )               +
      get(paste0("bedadpromedio_", alt))   * edad_promedio +#(n_trabajadores >= 2 )               +
      get(paste0("bq2_", alt)) * (quintil == 2) +
      get(paste0("bq3_", alt)) * (quintil == 3) +
      get(paste0("bq4_", alt)) * (quintil == 4) +
      get(paste0("bq5_", alt)) * (quintil == 5) +
      get(paste0("bnorte_", alt)) * (macrozona == "norte") +
      get(paste0("bcentro_", alt)) * (macrozona == "centro") +
      get(paste0("bsur_", alt)) * (macrozona == "sur")
  }

  ### Define settings for MNL models component
  fmnl_settings <- list(
    alternatives  = c("alimentos", "vestimenta", "cuentas", "hogar", "salud", "transporte","comunicaciones", "recreacion", "educacion", "restaurantes"),
    choiceShares  = list(alimentos = g.alimentos, vestimenta = g.vestimenta, cuentas = g.cuentas,
                         hogar = g.hogar, salud = g.salud, transporte = g.transporte, comunicaciones = g.comunicaciones,
                         recreacion = g.recreacion, educacion = g.educacion, restaurantes = g.restaurantes),
    utilities     = V)

  P[["model"]] <- apollo_fmnl(fmnl_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)
apollo_saveOutput(model)
apollo_loadModel("output/FMNL")
