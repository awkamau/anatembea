## ---------------------------------------------------------------------------
## Basic example: run an anatembea odin model using monthly EIR or betaa.
##
## This script shows two options:
##   1. Monthly EIR -> odin_model_stripped_matched.R
##   2. Monthly betaa -> MiP_odin_model_nodelay.R
##
## Both odin model files are loaded directly from the installed anatembea
## package. The monthly time/value vectors are passed directly to the odin
## model. The model handles interpolation internally and returns daily output.
## ---------------------------------------------------------------------------
#devtools::install_github("mrc-ide/anatembea")

library(anatembea)
library(odin)


## Age groups used by the anatembea transmission model.
age_vector <- c(
  0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5,
  5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80
)


## Compile one of the odin model files supplied with anatembea.
compile_anatembea_model <- function(model_file) {
  model_path <- system.file("odin", model_file, package = "anatembea")

  if (!nzchar(model_path)) {
    stop("Could not find ", model_file, " in the installed anatembea package.")
  }

  model_text <- readLines(model_path, warn = FALSE)

  ## do.call() is used because odin captures its model argument specially.
  do.call(odin::odin, list(model_text))
}


## Check a simple yearly input data frame.
##
## time is the model time for each yearly estimate. For example, if estimates
## are approximately yearly, this might be c(0, 365, 730, ...). Use the same
## time vector that accompanies the estimated EIR or betaa trajectory.
prepare_year_input <- function(year_data, value_column) {
  required_columns <- c("time", value_column)
  missing_columns <- setdiff(required_columns, names(year_data))

  if (length(missing_columns) > 0) {
    stop("Missing column(s): ", paste(missing_columns, collapse = ", "))
  }

  year_data <- year_data[, required_columns]
  year_data$time <- as.numeric(year_data$time)
  year_data[[value_column]] <- as.numeric(year_data[[value_column]])
  year_data <- year_data[order(year_data$time), ]

  if (anyNA(year_data) || any(!is.finite(year_data[[value_column]]))) {
    stop("Yearly times and values must not contain missing/non-finite values.")
  }
  if (anyDuplicated(year_data$time)) {
    stop("There must be only one value for each yearly time.")
  }

  year_data
}


## Create the model parameter list and equilibrium starting state.
##
## init_eir determines the equilibrium state at the beginning of the run.
## If log_init_EIR is available instead, use:
##   init_eir <- exp(log_init_EIR)
create_initial_state <- function(init_eir, model_arguments = list()) {
  parameter_arguments <- list(
    init_EIR = init_eir,
    comparison = "u5",
    lag_rates = 10
  )
  parameter_arguments <- c(parameter_arguments, model_arguments)

  model_parameters <- do.call(
    anatembea::model_param_list_create,
    parameter_arguments
  )

  equilibrium_state <- anatembea::equilibrium_init_create_stripped(
    age_vector = age_vector,
    het_brackets = 5,
    ft = 0.4,
    init_EIR = init_eir,
    model_param_list = model_parameters
  )

  list(
    model_parameters = model_parameters,
    equilibrium_state = equilibrium_state
  )
}


## Convert yearly EIR into daily mechanistic model output.
run_yearly_eir_daily <- function(yearly_eir, init_eir = NULL) {
  yearly_eir <- prepare_year_input(yearly_eir, "eir")

  if (is.null(init_eir)) {
    init_eir <- yearly_eir$eir[1]
  }

  eir_times <- yearly_eir$time
  eir_vals <- yearly_eir$eir

  initial_values <- create_initial_state(
    init_eir = init_eir,
    model_arguments = list(
      EIR_times = eir_times,
      EIR_vals = eir_vals
    )
  )

  generator <- compile_anatembea_model("odin_model_stripped_matched.R")

  ## Follow the same setup used in posterior_propagated_cases_averted.R.
  state_use <- initial_values$equilibrium_state[
    names(initial_values$equilibrium_state) %in% stats::coef(generator)$name
  ]
  state_use$lag_rates <- initial_values$model_parameters$lag_rates
  state_use$EIR_times <- eir_times
  state_use$EIR_vals <- eir_vals

  model <- generator$new(user = state_use, use_dde = TRUE)

  ## Only the requested output times are daily. The input trajectory remains
  ## yearly inside EIR_times and EIR_vals.
  tt <- seq(0, max(eir_times), by = 1)

  raw_output <- model$run(
    tt,
    step_max_n = 1e7,
    atol = 1e-5,
    rtol = 1e-5
  )
  output <- model$transform_variables(raw_output)

  ## EIR_eq is the internally interpolated EIR after age and biting-group
  ## scaling. Dividing out those scaling terms recovers the annual EIR
  ## trajectory supplied through EIR_times and EIR_vals.
  daily_input_eir <- as.numeric(output$EIR_eq[, 1, 1]) *
    state_use$DY /
    (state_use$rel_foi[1] * state_use$foi_age[1])

  data.frame(
    time = tt,
    input_eir_daily = daily_input_eir,
    ## EIR_out is calculated from the model's mosquito states:
    ## (av * Iv / omega) * DY. It is not the supplied EIR trajectory.
    mosquito_implied_eir = as.numeric(output$EIR_out),
    prev = as.numeric(output$prev),
    prev_all = as.numeric(output$prev_all),
    clin_inc = as.numeric(output$clin_inc),
    clinical_incidence_all_per_1000 = as.numeric(output$inc) * 1000,
    clinical_incidence_under5_per_1000 = as.numeric(output$incunder5) * 1000
  )
}


## Convert yearly betaa into daily mechanistic model output.
##
## The betaa model still needs init_eir to construct the equilibrium starting
## state. After initialization, the supplied betaa trajectory drives the model.
run_yearly_betaa_daily <- function(yearly_betaa, init_eir) {
  yearly_betaa <- prepare_year_input(yearly_betaa, "betaa")

  betaa_times <- yearly_betaa$time
  betaa_vals <- yearly_betaa$betaa

  initial_values <- create_initial_state(
    init_eir = init_eir,
    model_arguments = list(
      betaa_times = betaa_times,
      betaa_vals = betaa_vals
    )
  )

  generator <- compile_anatembea_model("MiP_odin_model_nodelay.R")

  state_use <- initial_values$equilibrium_state[
    names(initial_values$equilibrium_state) %in% stats::coef(generator)$name
  ]
  state_use$lag_rates <- initial_values$model_parameters$lag_rates
  state_use$betaa_times <- betaa_times
  state_use$betaa_vals <- betaa_vals

  model <- generator$new(user = state_use, use_dde = TRUE)

  tt <- seq(0, max(betaa_times), by = 1)

  raw_output <- model$run(
    tt,
    step_max_n = 1e7,
    atol = 1e-5,
    rtol = 1e-5
  )
  output <- model$transform_variables(raw_output)

  data.frame(
    time = tt,
    input_betaa = as.numeric(output$betaa_out),
    ## For the betaa-driven model, EIR_out is an output derived from the
    ## simulated mosquito states, so "mosquito_implied_eir" is appropriate.
    mosquito_implied_eir = as.numeric(output$EIR_out),
    clin_inc = as.numeric(output$clin_inc),
    clinical_incidence_all_per_1000 = as.numeric(output$inc) * 1000,
    clinical_incidence_under5_per_1000 = as.numeric(output$inc05) * 1000
  )
}


## ---------------------------------------------------------------------------
## Example input and use
## ---------------------------------------------------------------------------

## Replace these values with estimated yearly EIR values.
example_yearly_eir <- data.frame(
  time = c(0, 365, 730, 1095),
  eir = c(10, 20, 40, 30)
)

eir2 <- eir %>% mutate(time = (year -1990)*365) %>% dplyr::rename(eir = median) %>% dplyr::select(time, eir)

eir2 <- eir %>% mutate(time = (year -2007)*365) %>% dplyr::rename(eir = median) %>% dplyr::select(time, eir)


daily_from_eir <- run_yearly_eir_daily(
  yearly_eir = eir2,
  init_eir = eir2$eir[1]
)


head(daily_from_eir)
View(daily_from_eir)

yearly_est <- daily_from_eir %>%
  mutate(timestep =  as.integer(time/365) + 1990) %>%
  group_by(timestep) %>%
  summarise_at(vars(prev:clinical_incidence_under5_per_1000), sum, na.rm = TRUE)


plot(daily_from_eir$time, daily_from_eir$input_eir_daily, type = "l",
     xlab = "Time (days)", ylab = "Daily input EIR")
plot(daily_from_eir$time, daily_from_eir$prev, type = "l",
     xlab = "Time (days)", ylab = "Prevalence")
plot(daily_from_eir$time, daily_from_eir$clinical_incidence_all_per_1000, type = "l",
     xlab = "Time (days)", ylab = "Clinical Incidence (all ages per 1000)")



ggplot(data=daily_from_eir)+
  geom_line(aes(x=time,y=prev),color="#1F78B4") +
  geom_point(data=Kilifi_pr_det4,aes(x=((year-1990)*365),y=proportion),color='red') +
  geom_errorbar(data=Kilifi_pr_det4,aes(x=((year-1990)*365),ymin = lower, ymax = upper),
                width = 0.2, alpha = 0.5, linetype =  "dashed", color ="red") +
  scale_y_continuous(limits=c(0,0.7),expand=c(0,0))+
  labs(y='Prevalence under 5 years')+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggplot(data=daily_from_eir)+
  geom_line(aes(x=time,y=prev_all),color="#1F78B4")+
  geom_point(data=Kilifi_pr_det4,aes(x=((year-1990)*365),y=proportion),color='red') +
  geom_errorbar(data=Kilifi_pr_det4,aes(x=((year-1990)*365),ymin = lower, ymax = upper),
                width = 0.2, alpha = 0.5, linetype =  "dashed", color ="red") +
  scale_y_continuous(limits=c(0,0.7),expand=c(0,0))+
  labs(y='Prevalence all ages')+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggplot(data=daily_from_eir)+
  geom_line(aes(x=time,y=clinical_incidence_all_per_1000),color="brown") +
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))+
  labs(y="Clinical Incidence (all ages per 1000)")+
  theme_minimal()+
  theme(axis.title.x = element_blank())


ggplot(data=daily_from_eir)+
  geom_line(aes(x=time,y=clinical_incidence_under5_per_1000),color="brown")+
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))+
  labs(y="Clinical Incidence (under 5 years per 1000)")+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggplot(data=yearly_est)+
  geom_line(aes(x=timestep,y=clinical_incidence_under5_per_1000),color="brown")+
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))+
  labs(y="Clinical Incidence (under 5 years per 1000)")+
  theme_minimal()+
  theme(axis.title.x = element_blank())


ggplot(data=yearly_est)+
  geom_line(aes(x=timestep,y=clinical_incidence_all_per_1000),color="brown")+
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))+
  labs(y="Clinical Incidence (all ages per 1000)")+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggplot(data=yearly_est)+
  geom_line(aes(x=timestep,y=prev_all),color="#1F78B4")+
  geom_point(data=Kilifi_pr_det4,aes(x=year,y=proportion),color='red') +
  geom_errorbar(data=Kilifi_pr_det4,aes(x=year,ymin = lower, ymax = upper),
                width = 0.2, alpha = 0.5, linetype =  "dashed", color ="red") +
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))+
  labs(y="Prevalence under 5 years")+
  theme_minimal()+
  theme(axis.title.x = element_blank())


## Replace these values with estimated monthly betaa values.
example_monthly_betaa <- data.frame(
  time = c(0, 30, 60, 90),
  betaa = c(0.7, 0.9, 0.6, 0.8)
)

beta2 <- beta %>% mutate(time = (year -1990)*365) %>% dplyr::rename(betaa = median) %>% dplyr::select(time, betaa)

daily_from_betaa <- run_yearly_betaa_daily(
  yearly_betaa = beta2,
  init_eir = 12.314
)

head(daily_from_betaa)
plot(daily_from_betaa$time, daily_from_betaa$mosquito_implied_eir, type = "l",
     xlab = "Time (days)", ylab = "Mosquito-implied EIR")
plot(daily_from_betaa$time, daily_from_betaa$clinical_incidence_all_per_1000, type = "l",
     xlab = "Time (days)", ylab = "Clinical Incidence (all ages per 1000)")
plot(daily_from_betaa$time, daily_from_betaa$input_betaa, type = "l",
     xlab = "Time (days)", ylab = "Betaa")
