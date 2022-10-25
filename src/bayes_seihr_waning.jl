prob1 = ODEProblem(seir_ode_log!,
  zeros(11),
  (0.0, obstimes[end]),
  ones(10))

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
  dur_waning_non_centered_omicron ~ Normal()
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

  death_detection_rate_non_centered ~ Normal()

  ϕ_cases_non_centered ~ Exponential()
  ϕ_hospitalizations_non_centered ~ Normal()
  ϕ_death_non_centered ~ Exponential()
  ϕ_icu_non_centered ~ Normal()

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

  dur_icu_non_omicron = exp(dur_icu_non_centered_non_omicron * 0.1 + log(1))

  ω_non_omicron = 1/dur_icu_non_omicron

  dur_icu_omicron = exp(dur_icu_non_centered_omicron * 0.1 + log(0.31))
  ω_omicron = 1/dur_icu_omicron

  HICUR_non_omicron = logistic(HICUR_non_centered_non_omicron * 0.2 - 1.32)
  HICUR_omicron = logistic(HICUR_non_centered_omicron * 0.2 - 1.69)

  ICUDR_non_omicron = logistic(ICUDR_non_centered_non_omicron * 0.2 - 1.10)
  ICUDR_omicron = logistic(ICUDR_non_centered_omicron * 0.2 - 1.59)

  death_detection_rate = logistic(death_detection_rate_non_centered * 0.2 + 2.313635)
  ϕ_cases = ϕ_cases_non_centered^(-2)
  ϕ_hospitalizations = exp(ϕ_hospitalizations_non_centered * ϕ_hosp_sd + ϕ_hosp_mean)
  ϕ_death = ϕ_death_non_centered^(-2)
  ϕ_icu = exp(ϕ_icu_non_centered * ϕ_icu_sd + ϕ_icu_mean)

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
  D_non_omicron_init = 1
  D_omicron_init = 1
  R_init = 1
  C_non_omicron_init = I_non_omicron_init
  C_omicron_init = I_omicron_init
  S_init = popsize - (E_non_omicron_init + E_omicron_init + I_non_omicron_init + I_omicron_init + H_non_omicron_init + H_omicron_init + R_init + ICU_non_omicron_init + ICU_omicron_init + D_non_omicron_init + D_omicron_init)
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
  death_mean = death_detection_rate .* (sol_death_other + sol_death_omicron)


  for i in 1:l
    if (i < end_other_cases_time)
    data_est_other_cases[i] ~ NegativeBinomial2(max(other_cases_mean[i], 0.0), ϕ_cases)
    end
    data_est_omicron_cases[i] ~ NegativeBinomial2(max(omicron_cases_mean[i], 0.0), ϕ_cases)
    data_hospitalizations[i] ~ NegativeBinomial2(max(hospitalizations_mean[i], 0.0), ϕ_hospitalizations)
    data_icu[i] ~ NegativeBinomial2(max(icu_mean[i], 0.0), ϕ_icu)
    data_est_death[i] ~ NegativeBinomial2(max(death_mean[i], 0.0), ϕ_death)

  end

  β_t_non_omicron = vcat(β_init_non_omicron, β_t_values_no_init_non_omicron)
  β_t_omicron = vcat(β_init_omicron, β_t_values_no_init_omicron)
  R₀_t_non_omicron = β_t_non_omicron / ν_non_omicron
  R₀_t_omicron =  β_t_omicron / ν_omicron
  S_both = sol_reg_scale_array[1, :]
  S_omicron_only = sol_reg_scale_array[2, :]
  R₀_t = max.(β_t_omicron / ν_omicron, β_t_non_omicron / ν_non_omicron .* S_both / popsize)
  Rₜ_t = R₀_t .* (S_both + S_omicron_only) / popsize
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
    death_detection_rate = death_detection_rate,
    ϕ_death = ϕ_death,
    ϕ_cases = ϕ_cases,
    ϕ_hospitalizations = ϕ_hospitalizations,
    ϕ_icu = ϕ_icu,
    β_t_non_omicron = β_t_non_omicron,
    β_t_omicron = β_t_omicron,
    R₀_t_non_omicron = R₀_t_non_omicron,
    R₀_t_omicron = R₀_t_omicron,
    R₀_t = R₀_t,
    Rₜ_t = Rₜ_t,
    S_both = S_both,
    S_omicron_only = S_omicron_only,
    E_non_omicron = sol_reg_scale_array[3, :],
    E_omicron = sol_reg_scale_array[4, :],
    I_non_omicron = sol_reg_scale_array[5, :],
    I_omicron = sol_reg_scale_array[6, :],
    H_non_omicron = sol_reg_scale_array[7, :],
    H_omicron = sol_reg_scale_array[8, :],
    R = sol_reg_scale_array[9, :],
    C_non_omicron = sol_reg_scale_array[10, :],
    C_omicron = sol_reg_scale_array[11, :],
    sol_new_cases_other = sol_new_cases_other,
    sol_new_cases_omicron = sol_new_cases_omicron,
    sol_hospitalizations = sol_hospitalizations,
    other_cases_mean = other_cases_mean,
    omicron_cases_mean = omicron_cases_mean,
    hospitalizations_mean = hospitalizations_mean,
    ICU_non_omicron = sol_reg_scale_array[12, :],
    ICU_omicron = sol_reg_scale_array[13, :] ,
    D_non_omicron = sol_reg_scale_array[14, :],
    D_omicron = sol_reg_scale_array[15, :],
    icu_mean = icu_mean,
    death_mean = death_mean
  )
end;