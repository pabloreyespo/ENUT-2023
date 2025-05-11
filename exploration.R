rm(list = ls())  # .rs.restartR()
library(chilemapas)
set.seed(42)

source('utils.R')

getPackages <- function() {
  pkgs <- c("dplyr", "forcats", "geomtextpath", "ggplot2", "ggpubr", "glue", "gridExtra", "haven", "Hmisc", "lubridate", "purrr",
            "RColorBrewer", "readr", "reshape2", "stringr", "tibble", "tidyr", "see")
  invisible(lapply(pkgs, library, character.only=TRUE))
}

getPackages()
zval <- 3

MuchosColores <- c("#001219","#005F73","#0A9396","#94D2BD","#E9D8A6", "#EE9B00","#CA6702","#BB3E03","#AE2012","#9B2226")
#### --------------------------------------------------------------------------------------------------------------------
ipc = 0.362
get_data(free_activities = list(
            t_paid_work             = c("t_to"),
            t_leisure_socialization = c("t_vsyo_csar"),
            t_leisure_hobbies       = c("t_vsyo_aa") ,
            t_leisure_reading       = c("t_mcm_leer"),
            t_leisure_audio         = c("t_mcm_audio"),
            t_leisure_video         = c("t_mcm_video"),
            t_leisure_computer         = c("t_mcm_computador"),
            t_meals = c("t_cpag_comer"),
            t_commute = c("t_tt1"),
            t_commute_new = c("t_tt2"),
            t_sleep = c("t_cpag_dormir")
            ))
model_data <- model_data %>%  #%>% filter(ec > 0, ec <1)
  mutate(
    w  = w / (1 + ipc),
    Ec = Ec / (1 + ipc),
    I  = I / (1 + ipc),
    ec = ec / (1 + ipc)
  )
dbs = model_data
data_stack <- dbs %>%
  dplyr::select(id_persona, edad_anios, quintil, female, all_of(expenditures), all_of(times)) %>%
  melt(measure.vars= 8:18, variable.name = "tipo_tiempo") %>%
  mutate(value = as.numeric(value)) %>% rename(tiempo = value) %>%
  mutate(tipo_tiempo = case_when(
    tipo_tiempo == "t_paid_work"                ~ "Paid work",
    tipo_tiempo == "t_leisure_socialization" ~ "Socializing",
    tipo_tiempo == "t_leisure_hobbies"  ~ "Hobbies and sports",
    tipo_tiempo == "t_leisure_reading"        ~ "Reading Media",
    tipo_tiempo == "t_leisure_audio"        ~ "Listening Media",
    tipo_tiempo == "t_leisure_video"        ~ "Watching Media",
    tipo_tiempo == "t_leisure_computer"     ~ "Computer, Phone",
    tipo_tiempo == "t_meals" ~ "Meals",
    tipo_tiempo == "t_sleep" ~ "Sleep",
    tipo_tiempo == "t_commute" ~ "Commute",
    tipo_tiempo == "t_commute_new" ~ "Commute (new)",
    tipo_tiempo == "Tc" ~ "Other",
    T ~ "-"))

data_stack <- data_stack %>% mutate(gender_code = case_when(female == 1 ~ "Female", T ~ "Male"))

data_summary <- data_stack %>%
  mutate(participates = case_when(tiempo > 0 ~ 100, T ~ 0)) %>%
  group_by(tipo_tiempo, gender_code) %>%
  summarise(t_medio = mean(tiempo[tiempo > 0], na.rm = T), t_std = sd(tiempo[tiempo > 0], na.rm = T), participacion = mean(participates, na.rm = T)) %>%
  pivot_wider(names_from = "gender_code", values_from = c(participacion, t_medio, t_std), names_glue = "{gender_code}_{.value}")

data_summary_full <- data_stack %>%
  mutate(participates = case_when(tiempo > 0 ~ 100, T ~ 0)) %>%
  group_by(tipo_tiempo, gender_code) %>%
  summarise(t_medio = mean(tiempo, na.rm = T), t_std = sd(tiempo, na.rm = T), participacion = mean(participates, na.rm = T)) %>%
  pivot_wider(names_from = "gender_code", values_from = c(participacion, t_medio, t_std), names_glue = "{gender_code}_{.value}")


hist_tiempos <- data_stack %>%
  ggplot( aes(x=tiempo, fill = factor(tipo_tiempo, levels =c("Paid work", "Socializing", "Hobbies and sports", "Reading Media", "Listening Media", "Watching Media", "Computer, Phone", "Sleep", "Meals", "Commute", "Commute (new)", "Other")))) +
  geom_histogram(alpha=0.9, position = 'identity', color = "black") +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", n = 11, type = c('continuous'))) +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1", n = 11, type = c('continuous'))) +
  labs(x="Time allocated",y= "Frequency",fill = "Activities", caption = "ENUT 2022") +
  facet_wrap(. ~ factor(tipo_tiempo, levels =c("Paid work", "Socializing", "Hobbies and sports","Reading Media", "Listening Media", "Watching Media", "Computer, Phone", "Sleep", "Meals", "Commute", "Commute (new)", "Other")),
             scales = "free", nrow = 3) +
  theme_tesis
hist_tiempos

ggsave("plots/Time use histogram horizontal.png",
       plot = hist_tiempos, height = 6, width = 12)

hist_tiempos <- data_stack %>%
  filter(tipo_tiempo %in%  c("Socializing", "Hobbies and sports", "Reading Media", "Listening Media", "Watching Media", "Computer, Phone") ) %>%
  ggplot( aes(x=tiempo, fill = factor(tipo_tiempo, levels =c("Socializing", "Hobbies and sports", "Reading Media", "Listening Media", "Watching Media", "Computer, Phone")))) +
  geom_histogram(alpha=0.9, position = 'identity', color = "black") +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", n = 6, type = c('continuous'))) +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1", n = 6, type = c('continuous'))) +
  labs(x="Time allocated",y= "Frequency",fill = "Activities", caption = "ENUT 2015") +
  facet_wrap(. ~ factor(tipo_tiempo, levels =c("Socializing", "Hobbies and sports","Reading Media", "Listening Media", "Watching Media", "Computer, Phone")),
             scales = "free", nrow = 2) +
  theme_tesis
hist_tiempos

ggsave("plots/Time use leisure histogram horizontal.png",
       plot = hist_tiempos, height = 6, width = 12)

hist_wage <- dbs %>%
  ggplot( aes(x = w, fill = "Hourly wage")) +
  geom_histogram(alpha=0.8, position = 'identity', bins = 100) +
  scale_color_manual(values = MuchosColores) +
  scale_fill_manual(values = MuchosColores) +
  labs(x="1000 CLP",y= "Frequency", fill = "Variable" ,caption = "ENUT 2015") +
  geom_textvline(label = round(mean(dbs$w), 2), xintercept= mean(dbs$w), color = "red") +
  geom_textvline(label = round(mean(dbs$w) +  zval*sd(dbs$w), 2), xintercept= mean(dbs$w) +  zval*sd(dbs$w), color = "red", linetype="dotted") +
  geom_textvline(label = round(mean(dbs$w) -  zval*sd(dbs$w), 2), xintercept= mean(dbs$w)-  zval*sd(dbs$w), color = "red", linetype="dotted") +
  scale_x_continuous(limits = c( 0,30 )) +
  theme_tesis
hist_wage

ggsave("plots/Histogram wages.png",
       plot = hist_wage, height = 7, width = 7)
aux <- dbs %>% filter(I > 0)
hist_fixed_income <- aux %>%
  ggplot( aes(x = I, fill = "Fixed Income")) +
  geom_histogram(alpha=0.8, position = 'identity', bins = 100) +
  scale_color_manual(values = MuchosColores) +
  scale_fill_manual(values = MuchosColores) +
  labs(x="1000 CLP",y= "Frequency", fill = "Variable" ,caption = "ENUT 2015") +
  geom_text(label = paste0("n=", nrow(aux)), y = 50, x = 50, color = 'red') +
  geom_textvline(label = round(mean(aux$I), 2), xintercept= mean(aux$I), color = "red") +
  geom_textvline(label = round(mean(aux$I) +  zval*sd(aux$I), 2), xintercept= mean(aux$I) +  zval*sd(aux$I), color = "red", linetype="dotted") +
  geom_textvline(label = round(mean(aux$I) -  zval*sd(aux$I), 2), xintercept= mean(aux$I)-  zval*sd(aux$I), color = "red", linetype="dotted") +
  scale_x_continuous(limits = c( -10,200 )) +
  theme_tesis
hist_fixed_income

ggsave("plots/Histogram Fixed Income.png",
       plot = hist_fixed_income, height = 7, width = 7)

hist_available_time <- dbs %>%
  mutate(ta = 168 - Tc)  %>%
  ggplot( aes(x = ta, fill = "Available Time")) +
  geom_histogram(alpha=0.8, position = 'identity', bins = 100) +
  scale_color_manual(values = MuchosColores) +
  scale_fill_manual(values = MuchosColores) +
  labs(x="1000 CLP",y= "Frequency", fill = "Variable" ,caption = "ENUT 2015") +
  geom_textvline(label = round(mean(dbs$ta - dbs$Tc), 2), xintercept= mean(dbs$w- dbs$Tc), color = "red") +
  geom_textvline(label = round(mean(dbs$ta - dbs$Tc) +  zval*sd(dbs$ta - dbs$Tc), 2), xintercept= mean(dbs$ta - dbs$Tc) +  zval*sd(dbs$ta - dbs$Tc), color = "red", linetype="dotted") +
  geom_textvline(label = round(mean(dbs$ta - dbs$Tc) -  zval*sd(dbs$ta - dbs$Tc), 2), xintercept= mean(dbs$ta - dbs$Tc)-  zval*sd(dbs$ta - dbs$Tc), color = "red", linetype="dotted") +
  #scale_x_continuous(limits = c( -10,200 )) +
  theme_tesis
hist_available_time

ggsave("plots/Histogram Available Time.png",
       plot = hist_available_time, height = 7, width = 7)

hist_committed <- dbs %>%
  ggplot( aes(x = ec, fill = "Normalized committed expenses")) +
  geom_histogram(alpha=0.8, position = 'identity', bins = 100) +
  scale_color_manual(values = MuchosColores) +
  scale_fill_manual(values = MuchosColores) +
  geom_textvline(label = "ec=0", xintercept= 0, color = "black") +
  geom_textvline(label = round(mean(dbs$ec), 2), xintercept= mean(dbs$ec), color = "red") +
  geom_textvline(label = round(mean(dbs$ec) +  zval*sd(dbs$ec), 2), xintercept= mean(dbs$ec)+  zval*sd(dbs$ec), color = "red", linetype="dotted") +
  geom_textvline(label = round(mean(dbs$ec) -  zval*sd(dbs$ec), 2), xintercept= mean(dbs$ec)-  zval*sd(dbs$ec), color = "red", linetype="dotted") +
  labs(x="",y= "Frequency", fill = "Variable" ,caption = "ENUT 2015") +
  theme_tesis
hist_committed

ggsave("plots/Histogram committed expenses (ec chico).png",
       plot = hist_committed, height = 5, width = 5)

hist_committed <- dbs %>%
  ggplot( aes(x = Ec, fill = "Committed expenses")) +
  geom_histogram(alpha=0.8, position = 'identity', bins = 100) +
  scale_color_manual(values = MuchosColores) +
  scale_fill_manual(values = MuchosColores) +
  geom_textvline(label = "ec=0", xintercept= 0, color = "black") +
  geom_textvline(label = round(mean(dbs$Ec), 2), xintercept= mean(dbs$Ec), color = "red") +
  geom_textvline(label = round(mean(dbs$Ec) +  zval*sd(dbs$Ec), 2), xintercept= mean(dbs$Ec) +  zval*sd(dbs$Ec), color = "red", linetype="dotted") +
  geom_textvline(label = round(mean(dbs$Ec) -  zval*sd(dbs$Ec), 2), xintercept= mean(dbs$Ec)-  zval*sd(dbs$Ec), color = "red", linetype="dotted") +
  labs(x="",y= "Frequency", fill = "Variable" ,caption = "ENUT 2015") +
  theme_tesis
hist_committed

ggsave("plots/Histogram committed expenses.png",
       plot = hist_committed, height = 5, width = 5)

######################################################################################################################
########                                         CORRELACIONES                                                 #######
######################################################################################################################
#### INCOME VARIABLES
get_data()
original   <- c("edad_anios","mayor45","menor25", "underage_in_household", "only_worker",
                  "university", "n_menores6", "n_personas", "n_menores", "vive_pareja",
                  "fuentes_externas",
                  "w", "I", "ing_trab", "ing_personal", "ingreso_hogar", "metropolitana", 'norte', 'centro', 'sur')
replacement <- c("age","age>=45","age<25","# underage >=1","only worker",
                   "education= university degree",
                   "# children < 6 years old", "# people in household", "# underage in household",
                   "lives with partner",
                   "receives external help",
                   "wage rate", "earnings", "work ernings", "personal income", "household income", "metropolitan", 'north', 'center', 'south')
model_data[,replacement] <- model_data[,original]

socioecon          <- c("age", "# children < 6 years old", "wage rate", "personal income", "household income")
socioecon_binary   <- c("female", "age>=45", "# underage >=1", "only worker", "lives with partner", "university", "metropolitan", 'north', 'center', 'south', "receives external help")
#socioecon_binary   <- c("female", "age>=45", "# underage in household>=1", "only worker", "university", "Metropolitan")

### EXPLORAR tiempos
testRes = corrplot::cor.mtest(model_data[,c(times, expenditures)], conf.level = 0.95)
corr_AE <- cor(model_data[,c(times, expenditures)])
plot_corr_AE <- ggcorrplot::ggcorrplot(corr_AE, outline.col = "white", type = "lower",
                                         colors = c("#d55e00","white","#0071b2"), lab = TRUE) +
  labs(title = "time-use-expenditure variables") +
  theme_tesis + theme(legend.position = "right")
plot_corr_AE


corr_binary <- psych::tetrachoric(model_data[,socioecon_binary])
plot_corr_binary <- ggcorrplot::ggcorrplot(corr_binary$rho, outline.col = "white", type = "lower",
                                         colors = c("#d55e00","white","#0071b2"), lab = TRUE) +
  labs(title = "Socioeconomic binary variables (*)") +
  theme_tesis + theme(legend.position = "right", plot.title = element_text(colour="#2D4263", size = 15))
plot_corr_binary

corr_todo <- cor(model_data[,c(times, expenditures, socioecon_binary)])
plot_corr_todo <- ggcorrplot::ggcorrplot(corr_todo, outline.col = "white", type = "lower",
                                         colors = c("#d55e00","white","#0071b2"), lab = TRUE) +
  labs(title = "time-use-expenditure against socioeconomic variables") +
  theme_tesis + theme(legend.position = "right", plot.title = element_text(colour="#2D4263", size = 15))
plot_corr_todo

plot_corr_todo <- ggcorrplot::ggcorrplot(corr_todo[socioecon_binary, c(times, expenditures)], outline.col = "white", type = "lower",
                                         colors = c("#d55e00","white","#0071b2"), lab = TRUE) +
  labs(title = "time-use-expenditure vs socioeconomic") +
  theme_tesis + theme(legend.position = "right", plot.title = element_text(colour="#2D4263", size = 15))
plot_corr_todo

ggsave(paste0("plots/correlations_binary.png"), plot_corr_binary, dpi= 400, width=10, height=10)
figure <- ggarrange(
  plot_corr_AE, plot_corr_binary, plot_corr_todo,
          labels = c("A", "B", "C"), widths = c(1,1,1), ncol = 3) +
  labs(title = "Correlations")

ggsave(paste0("plots/correlations.png"), figure, dpi= 400, width=16, height=6)

#temp <- model_data %>% group_by(gender_segment) %>% summarise(mean_tot_expenses = mean(total_expenses))
#temp[2,2] / temp[1,2]






get_data(free_activities = list(
            t_paid_work             = c("t_to"),
            t_leisure_hobbies       = c("t_vsyo_aa") ,
            t_leisure_reading       = c("t_mcm_leer"),
            t_leisure_audio         = c("t_mcm_audio"),
            t_leisure_video         = c("t_mcm_video"),
            t_leisure_computer         = c("t_mcm_computador"),
            t_leisure_socialization = c("t_vsyo_csar")
            ), baskets = T)
table(model_data$class)