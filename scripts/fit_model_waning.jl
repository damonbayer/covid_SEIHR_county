using DrWatson
using Revise
using JLD2
using FileIO
using CSV
using DataFrames
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
using ForwardDiff
using Optim
using Random
using LineSearches

using covid_SEIHR_county

county_id =
    if length(ARGS) == 0
        1
    else
      parse(Int64, ARGS[1])
    end

priors_only = county_id == 0

if priors_only
  county_id = 1
end

if Sys.isapple() | Sys.iswindows()
  results_dir = "results/"
else
  results_dir = "//dfs6/pub/bayerd/covid_SEIHR_county/results/"
end

mkpath(results_dir)
mkpath(joinpath(results_dir, "posterior_predictive"))
mkpath(joinpath(results_dir, "generated_quantities"))
mkpath(joinpath(results_dir, "posterior_samples"))

savename_dict = Dict(:county_id => county_id)

time_interval_in_days = 7

county = subset(CSV.read("data/county_id_key.csv", DataFrame), :id => ByRow(x -> x == county_id))[1, :county]

max_t = 20
n_forecast_times = 4

dat = subset(CSV.read("data/cases_hospitalizations_by_county.csv", DataFrame), :county => ByRow(x -> x == county))
max_index = searchsortedlast(dat[!, :time], max_t)

initialization_values = subset(CSV.read("data/initialization_values.csv", DataFrame), :county => ByRow(x -> x == county))
popsize = subset(CSV.read("data/county_pop.csv", DataFrame), :County => ByRow(x -> x == county))[!, :Population][1]

news_cases_omicron_initial = initialization_values[1, :est_omicorn_cases]
new_cases_other_initial = initialization_values[1, :est_other_cases]
hospitalizations_initial = initialization_values[1, :hospitalizations]
icu_initial = initialization_values[1, :icu]

obstimes = float(dat[:, :time])
obstimes_forecast = vcat(obstimes, obstimes[end] .+ float.(1:n_forecast_times))

param_change_times = obstimes[1:(end - 1)]
param_change_times_forecast = obstimes_forecast[1:(end - 1)]


data_est_other_cases = dat[:, :est_other_cases]
data_est_other_cases_forecast = vcat(data_est_other_cases, repeat([data_est_other_cases[end]], n_forecast_times))
missing_est_other_cases_forecast = repeat([missing], length(data_est_other_cases) + n_forecast_times)

data_est_omicron_cases = dat[:, :est_omicron_cases]
data_est_omicron_cases_forecast = vcat(data_est_omicron_cases, repeat([data_est_omicron_cases[end]], n_forecast_times))
missing_est_omicron_cases_forecast = repeat([missing], length(data_est_omicron_cases) + n_forecast_times)

data_tests = dat[:, :tests]
data_tests_forecast = vcat(data_tests, repeat([data_tests[end]], n_forecast_times))

data_est_other_tests = dat[:, :est_other_tests]
data_est_other_tests_forecast = vcat(data_est_other_tests, repeat([data_est_other_tests[end]], n_forecast_times))
end_other_cases_time = findfirst(iszero, data_est_other_cases)

if isnothing(end_other_cases_time)
  end_other_cases_time = length(data_est_other_cases) + 1
end

data_est_omicron_tests = dat[:, :est_omicron_tests]
data_est_omicron_tests_forecast = vcat(data_est_omicron_tests, repeat([data_est_omicron_tests[end]], n_forecast_times))

data_hospitalizations = dat[:, :hospitalizations]
data_hospitalizations_forecast = vcat(data_hospitalizations, repeat([data_hospitalizations[end]], n_forecast_times))
missing_hospitalizations_forecast = repeat([missing], length(data_hospitalizations) + n_forecast_times)

data_icu = dat[:, :icu]
data_icu_forecast = vcat(data_icu, repeat([data_icu[end]], n_forecast_times))
missing_icu_forecast = repeat([missing], length(data_icu) + n_forecast_times)

data_est_death = dat[:, :est_deaths]
data_est_death_forecast = vcat(data_est_death, repeat([data_est_death[end]], n_forecast_times))
missing_est_death_forecast = repeat([missing], length(data_est_death) + n_forecast_times)




other_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.1)
omicron_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.6)
other_rho_samples = other_case_detection_rate_samples / median(data_est_other_tests)
omicron_rho_samples = omicron_case_detection_rate_samples / median(data_est_omicron_tests)

other_rho_mean = mean(logit.(other_rho_samples))
other_rho_sd = std(logit.(other_rho_samples))

omicron_rho_mean = mean(logit.(omicron_rho_samples))
omicron_rho_sd = std(logit.(omicron_rho_samples))

function NegativeBinomial2(μ, ϕ)
  p = logistic(log(ϕ) - log(μ))

  if p <= zero(p)
      p = nextfloat(zero(p))
  end

  if p >= one(p)
      p = prevfloat(one(p))
  end

  r = ϕ

  if r <= zero(r)
      r = nextfloat(zero(r))
  end

  Distributions.NegativeBinomial(r, p)
end

# Higher ϕ ⟹ lower variance
# var(NegativeBinomial2(100, 0.00001))
# var(NegativeBinomial2(100, 100000))

function seir_ode_log!(du, u, p, t)
  (S_both, S_omicron_only, E_non_omicron, E_omicron, I_non_omicron, I_omicron, H_non_omicron, H_omicron, R, C_non_omicron, C_omicron, 
  ICU_non_omicron, 
  ICU_omicron,
  D_non_omicron,
  D_omicron) = exp.(u)
  (β_init_non_omicron, β_init_omicron, γ_non_omicron, γ_omicron, ν_non_omicron, ν_omicron, η_non_omicron, η_omicron, IHR_non_omicron, IHR_omicron, κ,
HICUR_non_omicron, HICUR_omicron, ω_non_omicron, ω_omicron, ICUDR_non_omicron, ICUDR_omicron) = p
  N = S_both + S_omicron_only + E_non_omicron + E_omicron + I_non_omicron + I_omicron + H_non_omicron + H_omicron + R + ICU_non_omicron + ICU_omicron + D_omicron + D_non_omicron

  # infection = β * I * S / N
  infection_omicron_only_to_omicron = β_init_omicron * I_omicron * S_omicron_only / N
  infection_both_to_omicron = β_init_omicron * I_omicron * S_both / N
  infection_both_to_non_omicron = β_init_non_omicron * I_non_omicron * S_both / N
  
  infection_omicron = infection_both_to_omicron + infection_omicron_only_to_omicron
  infection_non_omicron = infection_both_to_non_omicron

  # progression = γ * E
  progression_non_omicron = γ_non_omicron * E_non_omicron
  progression_omicron = γ_omicron * E_omicron

  # hospitalization = ν * IHR * I
  hospitalization_non_omicron = ν_non_omicron * IHR_non_omicron * I_non_omicron
  hospitalization_omicron = ν_omicron * IHR_omicron * I_omicron

  # non_hospitalized_recovery = ν * (1 - IHR) * I
  non_hospitalized_recovery_non_omicron = ν_non_omicron * (1 - IHR_non_omicron) * I_non_omicron
  non_hospitalized_recovery_omicron = ν_omicron * (1 - IHR_omicron) * I_omicron

  # hospitalized_recovery = η * H
  hospitalized_recovery_non_omicron = η_non_omicron * (1-HICUR_non_omicron) * H_non_omicron
  hospitalized_recovery_omicron = η_omicron * (1-HICUR_omicron) * H_omicron

  # waning immunity = \kappa * R, and I will assume it only goes back to S_omicron_only
  waning_immunity = κ * R

  # progression from H to ICU_non_omicron
  icu_non_omicron = η_non_omicron * HICUR_non_omicron * H_non_omicron
  icu_omicron = η_omicron * HICUR_omicron * H_omicron

  # progression from ICU to deaths 
  death_non_omicron = ω_non_omicron * ICUDR_non_omicron * ICU_non_omicron
  death_omicron = ω_omicron * ICUDR_omicron * ICU_omicron 

  # progression from ICU to recovery 
  icu_recovery_non_omicron = ω_non_omicron * (1-ICUDR_non_omicron) * ICU_non_omicron
  icu_recovery_omicron = ω_omicron * (1-ICUDR_omicron) * ICU_omicron



  @inbounds begin
    du[1] = -(infection_both_to_omicron + infection_both_to_non_omicron) / S_both # S_both
    du[2] = (waning_immunity - infection_omicron_only_to_omicron) / S_omicron_only # S_omicron_only
    du[3] = (infection_non_omicron - progression_non_omicron) / E_non_omicron # E_non_omicron
    du[4] = (infection_omicron - progression_omicron) / E_omicron # E_omicron
    du[5] = (progression_non_omicron - (hospitalization_non_omicron + non_hospitalized_recovery_non_omicron)) / I_non_omicron # I_non_omicron
    du[6] = (progression_omicron - (hospitalization_omicron + non_hospitalized_recovery_omicron)) / I_omicron # I_omicron
    du[7] = (hospitalization_non_omicron - hospitalized_recovery_non_omicron) / H_non_omicron # H_non_omicron
    du[8] = (hospitalization_omicron - hospitalized_recovery_omicron) / H_omicron # H_omicron
    du[9] = (non_hospitalized_recovery_non_omicron + hospitalized_recovery_non_omicron + non_hospitalized_recovery_omicron + hospitalized_recovery_omicron + icu_recovery_non_omicron + 
    icu_recovery_omicron - waning_immunity) / R  # R
    du[10] = progression_non_omicron / C_non_omicron # C_non_omicron
    du[11] = progression_omicron / C_omicron # C_omicron
    du[12] = (icu_non_omicron - icu_recovery_non_omicron - death_non_omicron)/ICU_non_omicron
    du[13] = (icu_omicron - icu_recovery_omicron - death_omicron)/ICU_omicron
    du[14] = (death_non_omicron)/D_non_omicron
    du[15] = (death_omicron)/D_omicron
  end
  nothing
end

prob1 = ODEProblem(seir_ode_log!,
    zeros(11),
    (0.0, obstimes[end]),
    ones(10)
    )

@model function bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations, 
  data_est_other_tests, data_est_omicron_tests, data_icu, data_est_death, obstimes, param_change_times, extra_ode_precision, end_other_cases_time)
  l = length(data_est_other_cases)

  # Priors
  R0_params_non_centered ~ MvNormal(l + 3, 1) # +3 for sigma, non_omicron_init, omicron_init
  prop_omicron_only_init_non_centered ~ Normal()
  dur_latent_non_centered_non_omicron ~ Normal()
  dur_infectious_non_centered_non_omicron ~ Normal()
  IHR_non_centered_non_omicron ~ Normal()
  dur_hospitalized_non_centered_non_omicron ~ Normal()
  dur_icu_non_centered_non_omicron ~ Normal()
  HICUR_non_centered_non_omicron ~ Normal()
  ICUDR_non_centered_non_omicron ~ Normal()
  
  E_init_non_centered_non_omicron ~ Normal()
  I_init_non_centered_non_omicron ~ Normal()
  case_detection_rate_non_centered_other ~ Normal()
  
  dur_latent_non_centered_omicron ~ Normal()
  dur_infectious_non_centered_omicron ~ Normal()
  IHR_non_centered_omicron ~ Normal()
  dur_hospitalized_non_centered_omicron ~ Normal()
  dur_icu_non_centered_omicron ~ Normal()
  HICUR_non_centered_omicron ~ Normal()
  ICUDR_non_centered_omicron ~ Normal()
  dur_waning_non_centered_omicron ~ Normal()

  E_init_non_centered_omicron ~ Normal()
  I_init_non_centered_omicron ~ Normal()

  case_detection_rate_non_centered_omicron ~ Normal()

  ϕ_cases_non_centered ~ Exponential()
  ϕ_hospitalizations_non_centered ~ Exponential()
  ϕ_death_non_centered ~ Exponential()

  # troubleshooting 
  # R0_params_non_centered = zeros(l+3) # +3 for sigma, non_omicron_init, omicron_init
  # prop_omicron_only_init_non_centered = 0
  # dur_latent_non_centered_non_omicron = 0
  # dur_infectious_non_centered_non_omicron = 0
  # IHR_non_centered_non_omicron = 0
  # dur_hospitalized_non_centered_non_omicron = 0 
  # dur_icu_non_centered_non_omicron = 0
  # HICUR_non_centered_non_omicron = 0 
  # ICUDR_non_centered_non_omicron = 0
  
  # E_init_non_centered_non_omicron = 0 
  # I_init_non_centered_non_omicron = 0 
  # case_detection_rate_non_centered_other = 0
  
  # dur_latent_non_centered_omicron = 0
  # dur_infectious_non_centered_omicron = 0
  # IHR_non_centered_omicron = 0
  # dur_hospitalized_non_centered_omicron = 0
  # dur_icu_non_centered_omicron = 0
  # HICUR_non_centered_omicron = 0
  # ICUDR_non_centered_omicron = 0
  # dur_waning_non_centered_omicron = 0

  # E_init_non_centered_omicron = 0
  # I_init_non_centered_omicron = 0

  # case_detection_rate_non_centered_omicron  = 0

  # ϕ_cases_non_centered = 1 
  # ϕ_hospitalizations_non_centered  = 1
  # ϕ_death_non_centered = 1


  # Transformations
  R₀_init_non_centered_non_omicron = R0_params_non_centered[1]
  R₀_init_non_centered_omicron = R0_params_non_centered[2]
  σ_R0_non_centered = R0_params_non_centered[3]
  log_R0_steps_non_centered = R0_params_non_centered[4:end]

  R₀_init_non_omicron = exp(R₀_init_non_centered_non_omicron * 0.2 + 1.65)
  R₀_init_omicron = exp(R₀_init_non_centered_omicron * 0.2 + 0.55)
  σ_R0 = exp(σ_R0_non_centered * 0.2 - 3)
  
  case_detection_rate_other = logistic.(case_detection_rate_non_centered_other * 0.2 .- 1.1)
  case_detection_rate_omicron = logistic.(case_detection_rate_non_centered_omicron * 0.2 .- 1.6)

  dur_latent_non_omicron = exp(dur_latent_non_centered_non_omicron * 0.25 - 1.27)
  dur_latent_omicron = exp(dur_latent_non_centered_omicron * 0.25 - 1.5)

  γ_non_omicron = 1 / dur_latent_non_omicron
  γ_omicron = 1 / dur_latent_omicron

  dur_infectious_non_omicron = exp(dur_infectious_non_centered_non_omicron * 0.25 + log(6 / 7))
  dur_infectious_omicron = exp(dur_infectious_non_centered_omicron * 0.25 + log(5 / 7))

  ν_non_omicron = 1 / dur_infectious_non_omicron
  ν_omicron = 1 / dur_infectious_omicron

  β_init_non_omicron = R₀_init_non_omicron * ν_non_omicron
  β_init_omicron = R₀_init_omicron * ν_omicron

  IHR_non_omicron = logistic(IHR_non_centered_non_omicron * 0.2 - 3.2)
  IHR_omicron = logistic(IHR_non_centered_omicron * 0.25 - 4.3)

  dur_hospitalized_non_omicron = exp(dur_hospitalized_non_centered_non_omicron * 0.1 - 0.36)
  dur_hospitalized_omicron = exp(dur_hospitalized_non_centered_omicron * 0.1 - 1.54)

  η_non_omicron = 1 / dur_hospitalized_non_omicron
  η_omicron = 1 / dur_hospitalized_omicron

  dur_waning_omicron = exp(dur_waning_non_centered_omicron * 0.1 + log(12))

  κ = 1/dur_waning_omicron
# need priors from here to ICUDR_omicron 
  dur_icu_non_omicron = exp(dur_icu_non_centered_non_omicron * 0.1 + log(1))

  ω_non_omicron = 1/dur_icu_non_omicron

  dur_icu_omicron = exp(dur_icu_non_centered_omicron * 0.1 + log(0.31))
  ω_omicron = 1/dur_icu_omicron

  HICUR_non_omicron = logistic(HICUR_non_centered_non_omicron * 0.2 - 1.32)
  HICUR_omicron = logistic(HICUR_non_centered_omicron * 0.2 - 1.69)

  ICUDR_non_omicron = logistic(ICUDR_non_centered_non_omicron * 0.2 - 1.10)
  ICUDR_omicron = logistic(ICUDR_non_centered_omicron * 0.2 - 1.59)

  ϕ_cases = ϕ_cases_non_centered^(-2)
  ϕ_hospitalizations = ϕ_hospitalizations_non_centered^(-2)
  ϕ_death = ϕ_death_non_centered^(-2)

  # Initial state
  prop_omicron_only_init = logistic(prop_omicron_only_init_non_centered * 0.15 + 1.2)
  E_non_omicron_init = E_init_non_centered_non_omicron * 0.05 + new_cases_other_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  E_omicron_init = E_init_non_centered_omicron * 0.05 + news_cases_omicron_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  I_non_omicron_init = I_init_non_centered_non_omicron * 0.05 + new_cases_other_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  I_omicron_init = I_init_non_centered_omicron * 0.05 + news_cases_omicron_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  H_non_omicron_init = hospitalizations_initial
  H_omicron_init = 1
  ICU_non_omicron_init = icu_initial
  ICU_omicron_init = 1
  R_init = 1
  D_non_omicron_init = 1
  D_omicron_init = 1
  C_non_omicron_init = I_non_omicron_init
  C_omicron_init = I_omicron_init
  S_init = popsize - (E_non_omicron_init + E_omicron_init + I_non_omicron_init + I_omicron_init + H_non_omicron_init + H_omicron_init + R_init)
  S_both_init = (1 - prop_omicron_only_init) * S_init
  S_omicron_only_init = prop_omicron_only_init * S_init
  u0 = [S_both_init, S_omicron_only_init, E_non_omicron_init, E_omicron_init, I_non_omicron_init, I_omicron_init, H_non_omicron_init, H_omicron_init, R_init, C_non_omicron_init, C_omicron_init,
  ICU_non_omicron_init, ICU_omicron_init, D_non_omicron_init, D_omicron_init]

  # Time-varying parameters
  β_t_values_no_init_non_omicron = exp.(log(R₀_init_non_omicron) .+ cumsum(vec(log_R0_steps_non_centered) * σ_R0)) * ν_non_omicron
  β_t_values_no_init_omicron = exp.(log(R₀_init_omicron) .+ cumsum(vec(log_R0_steps_non_centered) * σ_R0)) * ν_omicron
  prob = remake(prob1, u0 = log.(u0), p = [β_init_non_omicron, β_init_omicron, γ_non_omicron, γ_omicron, ν_non_omicron, ν_omicron, η_non_omicron, η_omicron, IHR_non_omicron, IHR_omicron, κ,
  HICUR_non_omicron, HICUR_omicron, ω_non_omicron, ω_omicron, ICUDR_non_omicron, ICUDR_omicron], tspan = (0.0, obstimes[end]))
  
  function param_affect_β!(integrator)
    ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
    integrator.p[1] = β_t_values_no_init_non_omicron[ind_t] # Replace β with a new value from β_t_values_no_init_non_omicron
    integrator.p[2] = β_t_values_no_init_omicron[ind_t] # Replace β with a new value from β_t_values_no_init_omicron
  end
  
  param_callback = PresetTimeCallback(param_change_times, param_affect_β!, save_positions = (false, false))

  # extra_ode_precision = false
  if extra_ode_precision
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false, abstol = 1e-11, reltol = 1e-8)
  else
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)
  end
  
  if sol.retcode != :Success
    Turing.@addlogprob! -Inf
    return
  end
  sol_reg_scale_array = exp.(Array(sol))

  sol_new_cases_other = sol_reg_scale_array[10, 2:end] - sol_reg_scale_array[10, 1:(end-1)]
  sol_new_cases_omicron = sol_reg_scale_array[11, 2:end] - sol_reg_scale_array[11, 1:(end-1)]
  sol_hospitalizations = sol_reg_scale_array[7, 2:end] + sol_reg_scale_array[8, 2:end]
  sol_icu = sol_reg_scale_array[12, 2:end] + sol_reg_scale_array[13, 2:end]
  sol_death_other = sol_reg_scale_array[14, 2:end] - sol_reg_scale_array[14, 1:(end - 1)]
  sol_death_omicron = sol_reg_scale_array[15, 2:end] - sol_reg_scale_array[15, 1:(end - 1)]

  other_cases_mean = sol_new_cases_other .* case_detection_rate_other
  omicron_cases_mean = sol_new_cases_omicron .* case_detection_rate_omicron
  hospitalizations_mean = sol_hospitalizations
  icu_mean = sol_icu 
  death_mean = sol_death_other + sol_death_omicron 

  for i in 1:l
    if (i < end_other_cases_time)
      data_est_other_cases[i] ~ NegativeBinomial2(max(other_cases_mean[i], 0.0), ϕ_cases)
    end 
    data_est_omicron_cases[i] ~ NegativeBinomial2(max(omicron_cases_mean[i], 0.0), ϕ_cases)
    data_hospitalizations[i] ~ Poisson(max(hospitalizations_mean[i], 0.0))
    data_icu[i] ~ Poisson(max(icu_mean[i], 0.0))
    data_est_death[i] ~ NegativeBinomial2(max(death_mean[i], 0.0), ϕ_death)
  end
  return (
    σ_R0 = σ_R0,
    case_detection_rate_other = case_detection_rate_other,
    case_detection_rate_omicron = case_detection_rate_omicron,
    dur_latent_non_omicron_days = dur_latent_non_omicron * 7,
    dur_latent_omicron_days = dur_latent_omicron * 7,
    dur_infectious_non_omicron_days = dur_infectious_non_omicron * 7,
    dur_infectious_omicron_days = dur_infectious_omicron * 7,
    β_init_non_omicron = R₀_init_non_omicron * ν_non_omicron,
    β_init_omicron = R₀_init_omicron * ν_omicron,
    IHR_non_omicron = IHR_non_omicron,
    IHR_omicron = IHR_omicron,
    dur_hospitalized_non_omicron_days = dur_hospitalized_non_omicron * 7,
    dur_hospitalized_omicron_days = dur_hospitalized_omicron * 7,
    dur_waning_omicron_days = dur_waning_omicron * 7,
    dur_icu_non_omicron_days = dur_icu_non_omicron * 7,
    dur_icu_omicron_days = dur_icu_omicron * 7,
    HICUR_non_omicron = HICUR_non_omicron,
    HICUR_omicron = HICUR_omicron,
    ICUDR_non_omicron = ICUDR_non_omicron,
    ICUDR_omicron = ICUDR_omicron,
    ϕ_cases = ϕ_cases,
    ϕ_hospitalizations = ϕ_hospitalizations,
    β_t_non_omicron = vcat(β_init_non_omicron, β_t_values_no_init_non_omicron),
    β_t_omicron = vcat(β_init_omicron, β_t_values_no_init_omicron),
    R₀_t_non_omicron = vcat(β_init_non_omicron, β_t_values_no_init_non_omicron) / ν_non_omicron,
    R₀_t_omicron = vcat(β_init_omicron, β_t_values_no_init_omicron) / ν_omicron,
    S_both = sol_reg_scale_array[1, :],
    S_omicron_only = sol_reg_scale_array[2, :],
    E_non_omicron = sol_reg_scale_array[3, :],
    E_omicron = sol_reg_scale_array[4, :],
    I_non_omicron = sol_reg_scale_array[5, :],
    I_omicron = sol_reg_scale_array[6, :],
    H_non_omicron = sol_reg_scale_array[7, :],
    H_omicron = sol_reg_scale_array[8, :],
    R = sol_reg_scale_array[9, :],
    C_non_omicron = sol_reg_scale_array[10, :],
    C_omicron = sol_reg_scale_array[11, :],
    ICU_non_omicron = sol_reg_scale_array[12, :],
    ICU_omicron = sol_reg_scale_array[13, :] ,
    D_non_omicron = sol_reg_scale_array[14, :],
    D_omicron = sol_reg_scale_array[15, :],
    sol_new_cases_other = sol_new_cases_other,
    sol_new_cases_omicron = sol_new_cases_omicron,
    sol_hospitalizations = sol_hospitalizations,
    other_cases_mean = other_cases_mean,
    omicron_cases_mean = omicron_cases_mean,
    hospitalizations_mean = hospitalizations_mean,
    icu_mean = icu_mean,
    death_mean = death_mean
  )
end;

n_samples = 100
n_chains = 4

my_model = bayes_seihr(
  data_est_other_cases,
  data_est_omicron_cases,
  data_hospitalizations,
  data_est_other_tests,
  data_est_omicron_tests,
  data_icu, 
  data_est_death,
  obstimes,
  param_change_times,
  false,
  end_other_cases_time)

my_model_forecast = bayes_seihr(
  data_est_other_cases_forecast,
  data_est_omicron_cases_forecast,
  data_hospitalizations_forecast,
  data_est_other_tests_forecast,
  data_est_omicron_tests_forecast,
  data_icu_forecast, 
  data_est_death_forecast,
  obstimes_forecast,
  param_change_times_forecast,
  true,
  end_other_cases_time)

my_model_forecast_missing = bayes_seihr(
  missing_est_other_cases_forecast,
  missing_est_omicron_cases_forecast,
  missing_hospitalizations_forecast,
  data_est_other_tests_forecast,
  data_est_omicron_tests_forecast,
  missing_icu_forecast, missing_est_death_forecast,
  obstimes_forecast,
  param_change_times_forecast,
  true,
  end_other_cases_time)

# Sample Posterior
# Should sample prior sometimes

if priors_only
  Random.seed!(county_id)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
  Random.seed!(county_id)
  prior_samples_forecast_randn = augment_chains_with_forecast_samples(Chains(prior_samples, :parameters), my_model, my_model_forecast, "randn")
  prior_indices_to_keep = .!isnothing.(generated_quantities(my_model_forecast, prior_samples_forecast_randn));
  prior_samples_forecast_randn = ChainsCustomIndex(prior_samples_forecast_randn, prior_indices_to_keep);
  
  wsave(joinpath(results_dir, "prior_samples.jld2"), @dict prior_samples)
  
  Random.seed!(county_id)
  prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_forecast_randn)
  CSV.write(joinpath(results_dir, "prior_predictive.csv"), DataFrame(prior_predictive_randn))
  
  Random.seed!(county_id)
  gq_randn = get_gq_chains(my_model_forecast, prior_samples_forecast_randn);
  CSV.write(joinpath(results_dir, "prior_generated_quantities.csv"), DataFrame(gq_randn))
end

MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]

alg = Gibbs(NUTS(-1, 0.8, 
:prop_omicron_only_init_non_centered, 
:dur_latent_non_centered_non_omicron, 
:dur_infectious_non_centered_non_omicron, 
:IHR_non_centered_non_omicron, 
:dur_hospitalized_non_centered_non_omicron, 
:dur_waning_non_centered_omicron, 
:E_init_non_centered_non_omicron, 
:I_init_non_centered_non_omicron, 
:case_detection_rate_non_centered_other, 
:dur_latent_non_centered_omicron, 
:dur_infectious_non_centered_omicron, 
:IHR_non_centered_omicron, 
:dur_hospitalized_non_centered_omicron, 
:E_init_non_centered_omicron, 
:I_init_non_centered_omicron, 
:case_detection_rate_non_centered_omicron, 
:ϕ_cases_non_centered, 
:ϕ_hospitalizations_non_centered,
:ϕ_death, 
:dur_icu_omicron, 
:dur_icu_non_omicron, 
:HICUR_non_omicron, 
:HICUR_omicron, 
:ICUD_non_omicron, 
:ICUD_omicron),
ESS(:R0_params_non_centered))
Random.seed!(county_id)

# testing with just using zeros, doesn't work
# test_init = zeros(51)
# chain = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, init_params = test_init)

posterior_samples = sample(my_model, alg, MCMCThreads(), n_samples, n_chains, init_params = repeat([MAP_init], n_chains) .* collect(range(0.92, stop = 0.98, length = n_chains)))

posterior_samples_forecast_randn = augment_chains_with_forecast_samples(Chains(posterior_samples, :parameters), my_model, my_model_forecast, "randn")

indices_to_keep = .!isnothing.(generated_quantities(my_model_forecast, posterior_samples_forecast_randn));

posterior_samples_forecast_randn = ChainsCustomIndex(posterior_samples_forecast_randn, indices_to_keep);
wsave(joinpath(results_dir, "posterior_samples", savename("posterior_samples", savename_dict, "jld2")), @dict posterior_samples)

Random.seed!(county_id)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_forecast_randn)
CSV.write(joinpath(results_dir, "posterior_predictive", savename("posterior_predictive", savename_dict, "csv")), DataFrame(predictive_randn))

Random.seed!(county_id)
gq_randn = get_gq_chains(my_model_forecast, posterior_samples_forecast_randn);
CSV.write(joinpath(results_dir, "generated_quantities", savename("generated_quantities", savename_dict, "csv")), DataFrame(gq_randn))