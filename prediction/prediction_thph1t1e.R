rm(list = ls())  # .rs.restartR()
source('utils.R')
source('apollo_jaradiaz.R')
model <- apollo_loadModel('prediction/thph1t1e')

ipc = 0.362
get_data()
model_data <- model_data %>% filter(ec>0)

model_data <- model_data %>%
  mutate(Metropolitana = metropolitana) %>%
  mutate(
    w = w / (1 + ipc),
    Ec = Ec / (1 + ipc),
    ec = Ec / (w *(ta-Tc)),
      # db 2015
    wbar = 3.119154,
    Ecbar = 73.30111,
    tabar = 82.58632,
    wsd = 4.028217,
    Ecsd = 70.37009,
    tasd = 17.26111)

apollo_initialise()
database <- model_data

apollo_beta <- model$estimate
apollo_inputs = apollo_validateInputs(
  apollo_control = model$apollo_control,
  apollo_lcPars = model$apollo_lcPars,
  apollo_beta  = apollo_beta,
  apollo_fixed = c(),
  database = database,
  silent       = T)

pred = apollo_prediction(
  model = model,
  apollo_probabilities = model$apollo_probabilities,
  apollo_inputs = apollo_inputs,
  #prediction_settings=list(modelComponent = "model")
)
database[,"Tw_pred"] = pred[,3]

z = 3
database %>%
  filter(
    w > wbar-z*wsd, w < wbar+z*wsd,
    Ec > Ecbar-z*Ecsd, Ec < Ecbar+z*Ecsd
  ) %>%
  ggplot(aes(x = Tw, y = Tw_pred)) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  geom_point() +
  geom_smooth() +
  lims(x = c(0,max(database$Tw)+1), y = c(0,max(database$Tw_pred)+1))
