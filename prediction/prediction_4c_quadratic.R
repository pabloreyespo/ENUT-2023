rm(list = ls())  # .rs.restartR()
source('utils.R')
# model <- apollo_loadModel('prediction/4C_quadratic_tesis')
model <- apollo_loadModel('prediction/4C-quadratic')

ipc = 0.362
get_data(
  free_activities = list(
    Tw = c("t_to"),
    Tfleisure = c('t_vsyo_csar', 't_vsyo_aa', 't_mcm_leer', 't_mcm_audio', 't_mcm_video', 't_mcm_computador')),
  free_expenditures = list(Ef1 = c("alimentos","recreacion","restaurantes","comunicaciones","vestimenta")))
tabla = model_data[,all_activities] %>% colMeans() %>% as.data.frame()

model_data <- model_data %>%
  mutate(Metropolitana = metropolitana) %>%
  mutate(
    w = w / (1 + ipc),
    Ec = Ec / (1 + ipc),
    # db tesis
    #wbar = 2.985597,
    #Ecbar = 58.25054,
    #tabar = 84.08031,
    #wsd = 3.713072,
    #Ecsd = 56.95729,
    #tasd = 17.78265,

    # db 2015
    wbar = 3.119154,
    Ecbar = 73.30111,
    tabar = 82.58632,
    wsd = 4.028217,
    Ecsd = 70.37009,
    tasd = 17.26111,

    # db 2022
    #wbar = 4.19846,
    #Ecbar = 58.81324,
    #tabar = 75.74767,
    #wsd = 8.708877,
    #Ecsd = 78.09443,
    #tasd = 15.72352,

    ta= 168 - Tc)

model_data$wbar[1]
model_data$Ecbar[1]
model_data$tabar[1]

model_data$wsd[1]
model_data$Ecsd[1]
model_data$tasd[1]

apollo_initialise()
nClass = 4
database <- model_data
apollo_beta <- model$estimate
apollo_inputs = apollo_validateInputs(
  apollo_control = model$apollo_control,
  apollo_lcPars = model$apollo_lcPars,
  apollo_beta  = apollo_beta,
  apollo_fixed = c(),
  database = database,
  silent       = T)

apollo_attach(apollo_beta, apollo_inputs)
#on.exit(apollo_detach(apollo_beta, apollo_inputs))
res <- model$apollo_lcPars(apollo_beta, apollo_inputs)
#pi_values <- res$pi_values

res = apollo_prediction(
  model = model,
  apollo_probabilities = model$apollo_probabilities,
  apollo_inputs = apollo_inputs,
  #prediction_settings=list(modelComponent = "model")
)
pi =  model$apollo_lcPars(apollo_beta, apollo_inputs)$pi_values
database[,"Tw_pred"] = res$model[,3]

database[,c("Tw","Tw_pred")]


pl <- database %>% filter(Tw_pred <= 168, Tw_pred >=0) %>%
  ggplot(aes(x = Tw, y = Tw_pred)) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  geom_point() +
  geom_smooth()

z = 3
pl <- database %>%
  filter(
    w > wbar-z*wsd, w < wbar+z*wsd,
    ta > tabar-z*tasd, ta < tabar+z*tasd,
    Ec > Ecbar-z*Ecsd, Ec < Ecbar+z*Ecsd
  ) %>%
  ggplot(aes(x = Tw, y = Tw_pred)) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  geom_point() +
  geom_smooth()

ggsave('prediction/4c_quadratic.png', pl, width = 6,height =6)