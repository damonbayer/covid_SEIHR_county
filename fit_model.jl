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


news_cases_initial = initialization_values[1, :est_cases]
hospitalizations_initial = initialization_values[1, :hospitalizations]

function seir_ode_log!(du, u, p, t)
    (S, E, I, H, R, C) = exp.(u)
    (β, γ, ν, η, IHR) = p
    N = S + E + I + H + R
  
    infection = β * I * S / N
    progression = γ * E
    hospitalization = ν * IHR * I
    non_hospitalized_recovery = ν * (1 - IHR) * I
    hospitalized_recovery = η * H
  
    @inbounds begin
        du[1] = -infection / S # S
        du[2] = (infection - progression) / E # E
        du[3] = (progression - (hospitalization + non_hospitalized_recovery)) / I # I
        du[4] = (hospitalization - hospitalized_recovery) / H # H
        du[5] = (non_hospitalized_recovery + hospitalized_recovery) / R # R
        du[6] = progression / C # Cumulative Progressions
    end
    nothing
  end

data_new_cases = dat[:, :est_cases]
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

@model bayes_seihr(data_new_cases, data_hospitalizations) = begin
  l = length(data_new_cases)
  ode_eval_times = collect(range(time_interval_in_days / 7, step = time_interval_in_days / 7, length = l))
  # Priors
  R₀_non_centered ~ Normal()
  dur_latent_non_centered ~ Normal()
  dur_infectious_non_centered ~ Normal()
  IHR_non_centered ~ Normal()
  dur_hospitalized_non_centered ~ Normal()
  ϕ_cases_non_centered ~ Exponential()
  ϕ_hospitalizations_non_centered ~ Exponential()
  case_detection_rate_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  I_init_non_centered ~ Normal()

  # R₀ = exp(R₀_non_centered * 0.2 + 1.7)
  R₀ = exp(R₀_non_centered * 0.2 + 1)

  # dur_latent = exp(dur_latent_non_centered * 0.154 - 1.27)
  dur_latent = exp(dur_latent_non_centered * 0.25 - 1.27)
  γ = 1 / dur_latent
  
  # dur_infectious = exp(dur_infectious_non_centered * 0.128 - 0.568)
  dur_infectious = exp(dur_infectious_non_centered * 0.25 - 0.568)
  ν = 1 / dur_infectious

  β = R₀ * ν
  
  IHR = logistic(IHR_non_centered * 0.4 - 4.5)
  
  # dur_hospitalized = exp(dur_hospitalized_non_centered * 0.101 - 0.342)
  dur_hospitalized = exp(dur_hospitalized_non_centered * 0.25 - 0.342)
  η = 1 / dur_hospitalized

  ϕ_cases = ϕ_cases_non_centered^(-2)
  ϕ_hospitalizations = ϕ_hospitalizations_non_centered^(-2)
  case_detection_rate = logistic(case_detection_rate_non_centered * 0.2 - 1.4)


  E_init = E_init_non_centered * 0.5 + news_cases_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  I_init = I_init_non_centered * 0.5 + news_cases_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  # E_init = news_cases_t_0 * 5 / 3
  # I_init = news_cases_t_0 * 10 / 3
  H_init = hospitalizations_initial
  R_init = 1
  S_init = popsize - (E_init + I_init + H_init + R_init)
  
  u0 = [S_init, E_init, I_init, H_init, R_init, I_init]
  p = [β, γ, ν, η, IHR]
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
  latent_cases = ode_sol_array[6,2:end] - ode_sol_array[6,1:(end - 1)]
  latent_hospitalizations = ode_sol_array[4,2:end]
  
  for i in 1:l   
    data_new_cases[i] ~ NegativeBinomial2(max(latent_cases[i] * case_detection_rate, 0), ϕ_cases)
    data_hospitalizations[i] ~ NegativeBinomial2(max(latent_hospitalizations[i], 0), ϕ_hospitalizations)
  end
  return(R₀ = R₀, 
  dur_latent_days = dur_latent * 7, 
  dur_infectious_days = dur_infectious * 7,
  IHR = IHR,
  dur_hospitalized_days = dur_hospitalized * 7,
  case_detection_rate = case_detection_rate)
end;

n_samples = 2000
n_chains = 4
my_model = bayes_seihr(data_new_cases, data_hospitalizations)

univariate_param_names = [:R₀, :dur_latent_days, :dur_infectious_days, :IHR, :dur_hospitalized_days, :case_detection_rate]

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
# gq_posterior = generated_quantities(my_model, posterior_samples)
# gq_prior_chains = Chains(permutedims(reinterpret(reshape, Float64, NamedTuple{Tuple(univariate_param_names)}.(gq_posterior)), [2, 1, 3]), univariate_param_names)
posterior_predictive = predict(bayes_seihr(
  Vector{Union{Missing, Int64}}(undef, length(data_new_cases) + 5),
  Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations) + 5)),
  posterior_samples)

CSV.write(string(results_dir, "posterior_samples_", county_id, ".csv"), DataFrame(posterior_samples))
CSV.write(string(results_dir, "posterior_predictive_samples_", county_id, ".csv"), DataFrame(posterior_predictive))