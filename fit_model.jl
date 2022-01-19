county_id = parse(Int64, ARGS[1])
using CSV
using DataFrames
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
results_dir = "//dfs6/pub/bayerd/covid_SEIHR_county/results/"
mkpath(results_dir)

time_interval_in_days = 3

county = subset(CSV.read("data/county_id_key.csv", DataFrame), :id => ByRow(x -> x == county_id))[1, :county]

function NegativeBinomial2(μ, ϕ)
  p = 1 / (1 + μ / ϕ)
  r = ϕ

  Distributions.NegativeBinomial(r, p)
end

# Higher ϕ ⟹ lower variance
# var(NegativeBinomial2(100, 0.00001))
# var(NegativeBinomial2(100, 100000))

dat = subset(CSV.read("data/cases_hospitalizations_by_county.csv", DataFrame), :county => ByRow(x -> x == county))
initialization_values = subset(CSV.read("data/initialization_values.csv", DataFrame), :county => ByRow(x -> x == county))
popsize = subset(CSV.read("data/county_pop.csv", DataFrame), :County => ByRow(x -> x == county))[!, :Population][1]


# news_cases_initial = initialization_values[1, :est_cases]
news_cases_omicron_initial = initialization_values[1, :est_omicorn_cases]
news_cases_other_initial = initialization_values[1, :est_other_cases]
hospitalizations_initial = initialization_values[1, :hospitalizations]

function seir_ode_log!(du, u, p, t)
    (S, Eₙ, Eₒ, Iₙ, Iₒ, Hₙ, Hₒ, R, Cₙ, Cₒ) = exp.(u)
    (βₙ, βₒ, γₙ, γₒ, νₙ, νₒ, ηₙ, ηₒ, IHRₙ, IHRₒ) = p
    N = S + Eₙ + Eₒ + Iₙ + Iₒ + Hₙ + Hₒ + R

    # infection = β * I * S / N
    infectionₙ = βₙ * Iₙ * S / N
    infectionₒ = βₒ * Iₒ * S / N
    # progression = γ * E
    progressionₙ = γₙ * Eₙ
    progressionₒ = γₒ * Eₒ
    # hospitalization = ν * IHR * I
    hospitalizationₙ = νₙ * IHRₙ * Iₙ
    hospitalizationₒ = νₒ * IHRₒ * Iₒ
    # non_hospitalized_recovery = ν * (1 - IHR) * I
    non_hospitalized_recoveryₙ = νₙ * (1 - IHRₙ) * Iₙ
    non_hospitalized_recoveryₒ = νₒ * (1 - IHRₒ) * Iₒ
    # hospitalized_recovery = η * H
    hospitalized_recoveryₙ = ηₙ * Hₙ
    hospitalized_recoveryₒ = ηₒ * Hₒ

    @inbounds begin
      du[1] = -(infectionₙ + infectionₒ) / S # S
      du[2] = (infectionₙ - progressionₙ) / Eₙ # Eₙ
      du[3] = (infectionₒ - progressionₒ) / Eₒ # Eₒ
      du[4] = (progressionₙ - (hospitalizationₙ + non_hospitalized_recoveryₙ)) / Iₙ # Iₙ
      du[5] = (progressionₒ - (hospitalizationₒ + non_hospitalized_recoveryₒ)) / Iₒ # Iₒ
      du[6] = (hospitalizationₙ - hospitalized_recoveryₙ) / Hₙ # Hₙ
      du[7] = (hospitalizationₒ - hospitalized_recoveryₒ) / Hₒ # Hₒ
      du[8] = (non_hospitalized_recoveryₙ + hospitalized_recoveryₙ + non_hospitalized_recoveryₒ + hospitalized_recoveryₒ) / R  # R
      du[9] = progressionₙ / Cₙ # Cₙ
      du[10] = progressionₒ / Cₒ # Cₒ
    end
    nothing
  end

data_est_other_cases = dat[:, :est_other_cases]
data_est_omicron_cases = dat[:, :est_omicron_cases]
data_hospitalizations = dat[:, :hospitalizations]

# R₀_non_centered = 0
# dur_latent_non_centered = 0
# dur_infectious_non_centered = 0
# IHR_non_centered = 0
# dur_hospitalized_non_centered = 0
# ϕ_cases_non_centered = 1
# ϕ_hospitalizations_non_centered = 1
# case_detection_rate_non_centered = 0
# E_init_non_centered = 0
# I_init_non_centered = 0

@model bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations) = begin
  l = length(data_est_other_cases)
  ode_eval_times = collect(range(time_interval_in_days / 7, step = time_interval_in_days / 7, length = l))
  
  # Priors
  R₀_non_centeredₙ ~ Normal()
  dur_latent_non_centeredₙ ~ Normal()
  dur_infectious_non_centeredₙ ~ Normal()
  IHR_non_centeredₙ ~ Normal()
  dur_hospitalized_non_centeredₙ ~ Normal()
  E_init_non_centeredₙ ~ Normal()
  I_init_non_centeredₙ ~ Normal()
  
  R₀_non_centeredₒ ~ Normal()
  dur_latent_non_centeredₒ ~ Normal()
  dur_infectious_non_centeredₒ ~ Normal()
  IHR_non_centeredₒ ~ Normal()
  dur_hospitalized_non_centeredₒ ~ Normal()
  E_init_non_centeredₒ ~ Normal()
  I_init_non_centeredₒ ~ Normal()

  ϕ_cases_non_centered ~ Exponential()
  ϕ_hospitalizations_non_centered ~ Exponential()
  case_detection_rate_non_centered ~ Normal()
  
  R₀ₙ = exp(R₀_non_centeredₙ * 0.2 + 0.25)
  R₀ₒ = exp(R₀_non_centeredₒ * 0.2 + 1.25)

  dur_latentₙ = exp(dur_latent_non_centeredₙ * 0.25 - 1.27)
  dur_latentₒ = exp(dur_latent_non_centeredₒ * 0.25 - 1.5)

  γₙ = 1 / dur_latentₙ  
  γₒ = 1 / dur_latentₒ

  dur_infectiousₙ = exp(dur_infectious_non_centeredₙ * 0.25 - 0.568)
  dur_infectiousₒ = exp(dur_infectious_non_centeredₒ * 0.25 - 0.85)
  
  νₙ = 1 / dur_infectiousₙ
  νₒ = 1 / dur_infectiousₒ

  βₙ = R₀ₙ * νₙ
  βₒ = R₀ₒ * νₒ

  IHRₙ = logistic(IHR_non_centeredₙ * 0.2 - 4.2)
  IHRₒ = logistic(IHR_non_centeredₒ * 0.25 - 5.3)

  dur_hospitalizedₙ = exp(dur_hospitalized_non_centeredₙ * 0.1 - 0.36)
  dur_hospitalizedₒ = exp(dur_hospitalized_non_centeredₒ * 0.1 - 1.54)
  
  ηₙ = 1 / dur_hospitalizedₙ
  ηₒ = 1 / dur_hospitalizedₒ

  ϕ_cases = ϕ_cases_non_centered^(-2)
  ϕ_hospitalizations = ϕ_hospitalizations_non_centered^(-2)
  case_detection_rate = logistic(case_detection_rate_non_centered * 0.2 - 1.4)

  # Eₙ_init = E_init_non_centeredₙ * 0.5 + news_cases_other_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  # Eₒ_init  = E_init_non_centeredₒ * 0.5 + news_cases_omicron_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  # Iₙ_init  = I_init_non_centeredₙ * 0.5 + news_cases_other_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  # Iₒ_init  = I_init_non_centeredₒ * 0.5 + news_cases_omicron_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  Eₙ_init = E_init_non_centeredₙ * 0.05 + news_cases_other_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  Eₒ_init  = E_init_non_centeredₒ * 0.05 + news_cases_omicron_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  Iₙ_init  = I_init_non_centeredₙ * 0.05 + news_cases_other_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  Iₒ_init  = I_init_non_centeredₒ * 0.05 + news_cases_omicron_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  Hₙ_init  = hospitalizations_initial - 1
  Hₒ_init  = 1
  R_init  = 1
  Cₙ_init  = Iₙ_init
  Cₒ_init  = Iₒ_init
  S_init  = popsize - (Eₙ_init + Eₒ_init + Iₙ_init + Iₒ_init + Hₙ_init + Hₒ_init + R_init)
  
  u0 = [S_init, Eₙ_init, Eₒ_init, Iₙ_init, Iₒ_init, Hₙ_init, Hₒ_init, R_init, Cₙ_init, Cₒ_init]
  p = [βₙ, βₒ, γₙ, γₒ, νₙ, νₒ, ηₙ, ηₒ, IHRₙ, IHRₒ]
  tspan = (0.0, ode_eval_times[end])
  prob = ODEProblem(seir_ode_log!,
          log.(u0),
          tspan,
          p)

  sol = solve(prob, Tsit5(), saveat=ode_eval_times, save_start = true)

  if sol.retcode != :Success # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    Turing.@addlogprob! -Inf
    return
  end

  ode_sol_array = exp.(Array(sol))
  latent_cases_other = ode_sol_array[9,2:end] - ode_sol_array[9,1:(end - 1)]
  latent_cases_omicron = ode_sol_array[10,2:end] - ode_sol_array[10,1:(end - 1)]
  latent_hospitalizations = ode_sol_array[6,2:end] + ode_sol_array[7,2:end]

  for i in 1:l
    data_est_other_cases[i] ~ NegativeBinomial2(max(latent_cases_other[i] * case_detection_rate, 0), ϕ_cases)
    data_est_omicron_cases[i] ~ NegativeBinomial2(max(latent_cases_omicron[i] * case_detection_rate, 0), ϕ_cases)
    data_hospitalizations[i] ~ NegativeBinomial2(max(latent_hospitalizations[i], 0), ϕ_hospitalizations)
  end
  return(R₀ₙ = R₀ₙ,
  R₀ₒ = R₀ₒ,
  dur_latentₙ_days = dur_latentₙ * 7,
  dur_latentₒ_days = dur_latentₒ * 7,
  dur_infectiousₙ_days = dur_infectiousₙ * 7,
  dur_infectiousₒ_days = dur_infectiousₒ * 7,
  IHRₙ = IHRₙ,
  IHRₒ = IHRₒ,
  dur_hospitalizedₙ_days = dur_hospitalizedₙ * 7,
  dur_hospitalizedₒ_days = dur_hospitalizedₒ * 7,
  case_detection_rate = case_detection_rate)
end;

n_samples = 2000
n_chains = 4
my_model = bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations)

univariate_param_names = [:R₀ₙ, :R₀ₒ, :dur_latentₙ_days, :dur_latentₒ_days, :dur_infectiousₙ_days, :dur_infectiousₒ_days, :IHRₙ, :IHRₒ, :dur_hospitalizedₙ_days, :dur_hospitalizedₒ_days, :case_detection_rate]

if county_id == 1
  prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
  gq_prior = generated_quantities(my_model, prior_samples)
  gq_prior_chains = Chains(permutedims(reinterpret(reshape, Float64, NamedTuple{Tuple(univariate_param_names)}.(gq_prior)), [2, 1, 3]), univariate_param_names)
  CSV.write(string(results_dir, "prior_gq_samples.csv"), DataFrame(gq_prior_chains))
end

# prior_predictive = predict(bayes_seihr(
#   Vector{Union{Missing, Int64}}(undef, length(data_new_cases)),
#   Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations))),
#   prior_samples)


#Sample Posterior
posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains)
gq_posterior = generated_quantities(my_model, posterior_samples)
gq_posterior_chains = Chains(permutedims(reinterpret(reshape, Float64, (gq_posterior)), [2,1,3]), collect(keys(gq_posterior[1,1,1])))

posterior_predictive = predict(bayes_seihr(Vector{Union{Missing, Int64}}(undef, length(data_est_other_cases) + 12),
  Vector{Union{Missing, Int64}}(undef, length(data_est_omicron_cases) + 12),
  Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations) + 12)),
  posterior_samples)

CSV.write(string(results_dir, "posterior_samples_", county_id, ".csv"), DataFrame(posterior_samples))
CSV.write(string(results_dir, "posterior_predictive_samples_", county_id, ".csv"), DataFrame(posterior_predictive))