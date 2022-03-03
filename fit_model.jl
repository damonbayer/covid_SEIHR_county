county_id =
    if length(ARGS) == 0
        29
    else
      parse(Int64, ARGS[1])
    end
using CSV
using DataFrames
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
using ForwardDiff
using Optim

if Sys.isapple()
  results_dir = "results/"
else
  results_dir = "//dfs6/pub/bayerd/covid_SEIHR_county/results/"
end

mkpath(results_dir)

time_interval_in_days = 7

county = subset(CSV.read("data/county_id_key.csv", DataFrame), :id => ByRow(x -> x == county_id))[1, :county]

function NegativeBinomial2(μ, ϕ)

  # p = 1 / (1 + μ / ϕ)
  # p = ϕ / (ϕ + μ)
  # p = exp(log(ϕ) - log(ϕ + μ))
  p = logistic(log(ϕ) - log(μ))

  if p <= zero(p)
  #   println("p <= zero(p)")
    # println("μ: ", big(ForwardDiff.value(μ)))
    # println("ϕ: ", big(ForwardDiff.value(ϕ)))
    # println("p: ", big(ForwardDiff.value(p)))
  #   error("p <= zero(p)")
    p = nextfloat(zero(p))
  end

  if p >= one(p)
  #   println("p >= one(p)")
    # println("μ: ", big(ForwardDiff.value(μ)))
    # println("ϕ: ", big(ForwardDiff.value(ϕ)))
    # println("p: ", big(ForwardDiff.value(p)))
  #   error("p >= one(p)")
    p = prevfloat(one(p))
  end
  
  r = ϕ

  if r <= zero(r)
    r = nextfloat(zero(r))
  end

  # println("r")
  # println(r)
  # println("p")
  # println(p)

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

data_est_other_cases = dat[:, :est_other_cases]
data_est_omicron_cases = dat[:, :est_omicron_cases]
data_tests = dat[:, :tests]
data_est_other_tests = dat[:, :est_other_tests]
data_est_omicron_tests = dat[:, :est_omicron_tests]
data_hospitalizations = dat[:, :hospitalizations]

case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.1)
rho_samples = case_detection_rate_samples / data_tests
rho_samples = case_detection_rate_samples / median(data_tests)
rho_mean = mean(logit.(rho_samples))
rho_sd = std(logit.(rho_samples))

other_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.1)
omicron_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.6)
other_rho_samples = other_case_detection_rate_samples / median(data_est_other_tests)
omicron_rho_samples = omicron_case_detection_rate_samples / median(data_est_omicron_tests)

other_rho_mean = mean(logit.(other_rho_samples))
other_rho_sd = std(logit.(other_rho_samples))

omicron_rho_mean = mean(logit.(omicron_rho_samples))
omicron_rho_sd = std(logit.(omicron_rho_samples))

function seir_ode_log!(du, u, p, t)
    (S_both, S_omicron_only, E_non_omicron, E_omicron, I_non_omicron, I_omicron, H_non_omicron, H_omicron, R, C_non_omicron, C_omicron) = exp.(u)
    (β_non_omicron, β_omicron, γ_non_omicron, γ_omicron, ν_non_omicron, ν_omicron, η_non_omicron, η_omicron, IHR_non_omicron, IHR_omicron) = p
    N = S_both + S_omicron_only + E_non_omicron + E_omicron + I_non_omicron + I_omicron + H_non_omicron + H_omicron + R

    # infection = β * I * S / N
    infection_omicron_only_to_omicron = β_omicron * I_omicron * S_omicron_only / N
    infection_both_to_omicron = β_omicron * I_omicron * S_both / N
    infection_both_to_non_omicron = β_non_omicron * I_non_omicron * S_both / N
    
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
    hospitalized_recovery_non_omicron = η_non_omicron * H_non_omicron
    hospitalized_recovery_omicron = η_omicron * H_omicron

    @inbounds begin
      du[1] = -(infection_both_to_omicron + infection_both_to_non_omicron) / S_both # S_both
      du[2] = -(infection_omicron_only_to_omicron) / S_omicron_only # S_omicron_only
      du[3] = (infection_non_omicron - progression_non_omicron) / E_non_omicron # E_non_omicron
      du[4] = (infection_omicron - progression_omicron) / E_omicron # E_omicron
      du[5] = (progression_non_omicron - (hospitalization_non_omicron + non_hospitalized_recovery_non_omicron)) / I_non_omicron # I_non_omicron
      du[6] = (progression_omicron - (hospitalization_omicron + non_hospitalized_recovery_omicron)) / I_omicron # I_omicron
      du[7] = (hospitalization_non_omicron - hospitalized_recovery_non_omicron) / H_non_omicron # H_non_omicron
      du[8] = (hospitalization_omicron - hospitalized_recovery_omicron) / H_omicron # H_omicron
      du[9] = (non_hospitalized_recovery_non_omicron + hospitalized_recovery_non_omicron + non_hospitalized_recovery_omicron + hospitalized_recovery_omicron) / R  # R
      du[10] = progression_non_omicron / C_non_omicron # C_non_omicron
      du[11] = progression_omicron / C_omicron # C_omicron
    end
    nothing
  end

l = length(data_est_other_cases)

# @model bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations, data_est_other_tests, data_est_omicron_tests) = begin
@model bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations, data_tests) = begin
  l = length(data_est_other_cases)
  ode_eval_times = collect(range(time_interval_in_days / 7, step = time_interval_in_days / 7, length = l))
  
  # Priors
  prop_omicron_only_init_non_centered ~ Normal()
  R₀_non_centered_non_omicron ~ Normal()
  dur_latent_non_centered_non_omicron ~ Normal()
  dur_infectious_non_centered_non_omicron ~ Normal()
  IHR_non_centered_non_omicron ~ Normal()
  dur_hospitalized_non_centered_non_omicron ~ Normal()
  E_init_non_centered_non_omicron ~ Normal()
  I_init_non_centered_non_omicron ~ Normal()
  # ρ_non_centered_non_omicron ~ Normal()
  case_detection_rate_non_centered_other ~ Normal()

  R₀_non_centered_omicron ~ Normal()
  dur_latent_non_centered_omicron ~ Normal()
  dur_infectious_non_centered_omicron ~ Normal()
  IHR_non_centered_omicron ~ Normal()
  dur_hospitalized_non_centered_omicron ~ Normal()
  E_init_non_centered_omicron ~ Normal()
  I_init_non_centered_omicron ~ Normal()
  # ρ_non_centered_omicron ~ Normal()
  case_detection_rate_non_centered_omicron ~ Normal()

  # ρ_non_centered ~ Normal()
  ϕ_cases_non_centered ~ Exponential()
  ϕ_hospitalizations_non_centered ~ Exponential()
  
  case_detection_rate_other = logistic.(case_detection_rate_non_centered_other * 0.2 .- 1.1)
  case_detection_rate_omicron = logistic.(case_detection_rate_non_centered_omicron * 0.2 .- 1.6)

  # Transformed Quantities
  prop_omicron_only_init = logistic(prop_omicron_only_init_non_centered * 0.15  + 1.2)

  # R₀_non_omicron = exp(R₀_non_centered_non_omicron * 0.2 + 0.25)
  # R₀_non_omicron = exp(R₀_non_centered_non_omicron * 0.2 + 1.5)
  R₀_non_omicron = exp(R₀_non_centered_non_omicron * 0.2 + 1.65)

  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 1.25)
  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 1.5)
  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 0)
  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 0.5) # This actually seems to work pretty well
  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 0.75) # This actually seems to work pretty well
  # R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 0.6) # This actually seems to work pretty well
  R₀_omicron = exp(R₀_non_centered_omicron * 0.2 + 0.55) # This actually seems to work pretty well

  dur_latent_non_omicron = exp(dur_latent_non_centered_non_omicron * 0.25 - 1.27)
  dur_latent_omicron = exp(dur_latent_non_centered_omicron * 0.25 - 1.5)

  γ_non_omicron = 1 / dur_latent_non_omicron  
  γ_omicron = 1 / dur_latent_omicron

  # dur_infectious_non_omicron = exp(dur_infectious_non_centered_non_omicron * 0.25 - 0.568)
  # dur_infectious_non_omicron = exp(dur_infectious_non_centered_non_omicron * 0.25 - 0.336)
  dur_infectious_non_omicron = exp(dur_infectious_non_centered_non_omicron * 0.25 + log(6/7))
  # dur_infectious_omicron = exp(dur_infectious_non_centered_omicron * 0.25 - 0.85)
  # dur_infectious_omicron = exp(dur_infectious_non_centered_omicron * 0.25 - 0.559)
  dur_infectious_omicron = exp(dur_infectious_non_centered_omicron * 0.25 + log(5/7))
  
  ν_non_omicron = 1 / dur_infectious_non_omicron
  ν_omicron = 1 / dur_infectious_omicron

  β_non_omicron = R₀_non_omicron * ν_non_omicron
  β_omicron = R₀_omicron * ν_omicron

  # IHR_non_omicron = logistic(IHR_non_centered_non_omicron * 0.2 - 4.2)
  IHR_non_omicron = logistic(IHR_non_centered_non_omicron * 0.2 - 3.2)
  # IHR_omicron = logistic(IHR_non_centered_omicron * 0.25 - 5.3)
  IHR_omicron = logistic(IHR_non_centered_omicron * 0.25 - 4.3)


  dur_hospitalized_non_omicron = exp(dur_hospitalized_non_centered_non_omicron * 0.1 - 0.36)
  dur_hospitalized_omicron = exp(dur_hospitalized_non_centered_omicron * 0.1 - 1.54)
  
  η_non_omicron = 1 / dur_hospitalized_non_omicron
  η_omicron = 1 / dur_hospitalized_omicron

  # ρ_non_omicron = exp(ρ_non_centered_non_omicron * other_rho_sd + other_rho_mean)
  # ρ_omicron = exp(ρ_non_centered_omicron *  omicron_rho_sd + omicron_rho_mean)
  # ρ = exp(ρ_non_centered * rho_sd + rho_mean)

  ϕ_cases = ϕ_cases_non_centered^(-2)
  ϕ_hospitalizations = ϕ_hospitalizations_non_centered^(-2)
  # case_detection_rate = logistic(case_detection_rate_non_centered * 0.2 - 1.4)

  E_non_omicron_init = E_init_non_centered_non_omicron * 0.05 + news_cases_other_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  E_omicron_init  = E_init_non_centered_omicron * 0.05 + news_cases_omicron_initial * 5 / 3 # new_cases_in_week_prior_to_model_start * (5/6, 20/6)
  I_non_omicron_init  = I_init_non_centered_non_omicron * 0.05 + news_cases_other_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  I_omicron_init  = I_init_non_centered_omicron * 0.05 + news_cases_omicron_initial * 10 / 3 # new_cases_in_week_prior_to_model_start * (5/3, 20/3)
  H_non_omicron_init  = hospitalizations_initial - 1
  H_omicron_init  = 1
  R_init  = 1
  C_non_omicron_init  = I_non_omicron_init
  C_omicron_init  = I_omicron_init
  S_init  = popsize - (E_non_omicron_init + E_omicron_init + I_non_omicron_init + I_omicron_init + H_non_omicron_init + H_omicron_init + R_init)
  S_both_init = (1 - prop_omicron_only_init) * S_init
  S_omicron_only_init = prop_omicron_only_init * S_init
  
  u0 = [S_both_init, S_omicron_only_init, E_non_omicron_init, E_omicron_init, I_non_omicron_init, I_omicron_init, H_non_omicron_init, H_omicron_init, R_init, C_non_omicron_init, C_omicron_init]
  p = [β_non_omicron, β_omicron, γ_non_omicron, γ_omicron, ν_non_omicron, ν_omicron, η_non_omicron, η_omicron, IHR_non_omicron, IHR_omicron]
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
  latent_cases_other = ode_sol_array[10,2:end] - ode_sol_array[10,1:(end - 1)]
  latent_cases_omicron = ode_sol_array[11,2:end] - ode_sol_array[11,1:(end - 1)]
  latent_hospitalizations = ode_sol_array[7,2:end] + ode_sol_array[8,2:end]
  
  # other_cases_mean = max.(latent_cases_other .* data_tests * ρ, repeat([0.0], l))
  # omicron_cases_mean = max.(latent_cases_omicron .* data_tests * ρ, repeat([0.0], l))

  # other_cases_mean = max.(latent_cases_other .* data_tests *ρ_non_omicron, repeat([0.0], l))
  # omicron_cases_mean = max.(latent_cases_omicron .* data_tests *ρ_omicron, repeat([0.0], l))

  other_cases_mean = max.(latent_cases_other .* case_detection_rate_other, repeat([0.0], l))
  omicron_cases_mean = max.(latent_cases_omicron .* case_detection_rate_omicron, repeat([0.0], l))

  if any(isnan.(other_cases_mean))
    println("WE GOT A NAN other_cases_mean")
    return
  end

  if any(isnan.(omicron_cases_mean))
    println("WE GOT A NAN omicron_cases_mean")
    return
  end

  for i in 1:l
    # data_est_other_cases[i] ~ NegativeBinomial2(latent_cases_other[i] * data_est_other_tests[i] * ρ_non_omicron, ϕ_cases)
    # data_est_omicron_cases[i] ~ NegativeBinomial2(latent_cases_omicron[i] * data_est_omicron_tests[i] * ρ_omicron, ϕ_cases)
    # println("data_est_other_cases")
    data_est_other_cases[i] ~ NegativeBinomial2(other_cases_mean[i], ϕ_cases)
    # println("data_est_omicron_cases") 
    data_est_omicron_cases[i] ~ NegativeBinomial2(omicron_cases_mean[i], ϕ_cases)
    # println("data_hospitalizations")
    data_hospitalizations[i] ~ NegativeBinomial2(max(latent_hospitalizations[i], 0.0), ϕ_hospitalizations)
  end
  return(prop_omicron_only_init = prop_omicron_only_init,
  R₀_non_omicron = R₀_non_omicron,
  R₀_omicron = R₀_omicron,
  dur_latent_non_omicron_days = dur_latent_non_omicron * 7,
  dur_latent_omicron_days = dur_latent_omicron * 7,
  dur_infectious_non_omicron_days = dur_infectious_non_omicron * 7,
  dur_infectious_omicron_days = dur_infectious_omicron * 7,
  IHR_non_omicron = IHR_non_omicron,
  IHR_omicron = IHR_omicron,
  dur_hospitalized_non_omicron_days = dur_hospitalized_non_omicron * 7,
  dur_hospitalized_omicron_days = dur_hospitalized_omicron * 7,
  ϕ_cases_non_centered = ϕ_cases_non_centered,
  ϕ_hospitalizations_non_centered = ϕ_hospitalizations_non_centered,
  # ρ_omicron = ρ_omicron,
  # ρ_non_omicron = ρ_non_omicron,
  # ρ = ρ,
  case_detection_rate_omicron = case_detection_rate_omicron,
  case_detection_rate_other = case_detection_rate_other,
  S_both = ode_sol_array[1,:],
  S_omicron_only = ode_sol_array[2,:],
  E_non_omicron = ode_sol_array[3,:],
  E_omicron = ode_sol_array[4,:],
  I_non_omicron = ode_sol_array[5,:],
  I_omicron = ode_sol_array[6,:],
  H_non_omicron = ode_sol_array[7,:],
  H_omicron = ode_sol_array[8,:],
  R = ode_sol_array[9,:],
  latent_cases_other = latent_cases_other,
  latent_cases_omicron = latent_cases_omicron,
  other_cases_mean = other_cases_mean,
  omicron_cases_mean = omicron_cases_mean)
end;

n_samples = 2000
n_chains = 4
my_model = bayes_seihr(data_est_other_cases, data_est_omicron_cases, data_hospitalizations, data_tests)

# univariate_param_names = [:prop_omicron_only_init, :R₀_non_omicron, :R₀_omicron, :dur_latent_non_omicron_days, :dur_latent_omicron_days, :dur_infectious_non_omicron_days, :dur_infectious_omicron_days, :IHR_non_omicron, :IHR_omicron, :dur_hospitalized_non_omicron_days, :dur_hospitalized_omicron_days, :ϕ_cases_non_centered, :ϕ_hospitalizations_non_centered, :ρ]
univariate_param_names = [:prop_omicron_only_init, :R₀_non_omicron, :R₀_omicron, :dur_latent_non_omicron_days, :dur_latent_omicron_days, :dur_infectious_non_omicron_days, :dur_infectious_omicron_days, :IHR_non_omicron, :IHR_omicron, :dur_hospitalized_non_omicron_days, :dur_hospitalized_omicron_days, :ϕ_cases_non_centered, :ϕ_hospitalizations_non_centered, :case_detection_rate_omicron, :case_detection_rate_other]
# univariate_param_names = [:prop_omicron_only_init, :R₀_non_omicron, :R₀_omicron, :dur_latent_non_omicron_days, :dur_latent_omicron_days, :dur_infectious_non_omicron_days, :dur_infectious_omicron_days, :IHR_non_omicron, :IHR_omicron, :dur_hospitalized_non_omicron_days, :dur_hospitalized_omicron_days, :ϕ_cases_non_centered, :ϕ_hospitalizations_non_centered, :ρ_omicron, :ρ_non_omicron]
time_varying_param_names = [:S_both, :S_omicron_only, :E_non_omicron, :E_omicron, :I_non_omicron, :I_omicron, :H_non_omicron, :H_omicron, :R, :other_cases_mean, :omicron_cases_mean, :latent_cases_other, :latent_cases_omicron]

function get_time_varying_param_chain(gq_chain, param_name)
  l = length(gq_chain[1][Symbol(param_name)]) - 1
  Chains(permutedims(reshape(hcat(getindex.(gq_chain, Symbol(param_name))...), l + 1, n_samples, n_chains), [2, 1, 3]),
  Symbol.(string.(string(param_name,"["), string.(0:l), "]")))
end

if county_id == 1
  Random.seed!(1)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
  gq_prior = generated_quantities(my_model, prior_samples);
  gq_prior_chains = hcat(vcat(Chains(permutedims(reinterpret(reshape, Float64, NamedTuple{Tuple(univariate_param_names)}.(gq_prior)), [2, 1, 3]), univariate_param_names),
  get_time_varying_param_chain.([gq_prior], time_varying_param_names))...)

  prior_predictive = predict(bayes_seihr(
    Vector{Union{Missing, Int64}}(undef, length(data_est_other_cases) + 12),
    Vector{Union{Missing, Int64}}(undef, length(data_est_omicron_cases) + 12),
    Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations) + 12),
    vcat(data_tests, repeat(data_tests[(end - 3):end], 3))),
    prior_samples)
  
  CSV.write(string(results_dir, "prior_samples.csv"), DataFrame(prior_samples))
  CSV.write(string(results_dir, "prior_gq_samples.csv"), DataFrame(gq_prior_chains))
  CSV.write(string(results_dir, "prior_predictive.csv"), DataFrame(prior_predictive))
end

#Sample Posterior
Random.seed!(county_id)
posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains)
CSV.write(string(results_dir, "posterior_samples_", county_id, ".csv"), DataFrame(posterior_samples))

gq_posterior = generated_quantities(my_model, posterior_samples)

gq_posterior_chains = hcat(vcat(Chains(permutedims(reinterpret(reshape, Float64, NamedTuple{Tuple(univariate_param_names)}.(gq_posterior)), [2, 1, 3]), univariate_param_names),
get_time_varying_param_chain.([gq_posterior], time_varying_param_names))...)

CSV.write(string(results_dir, "posterior_gq_samples_", county_id, ".csv"), DataFrame(gq_posterior_chains))

posterior_predictive = predict(bayes_seihr(
  Vector{Union{Missing, Int64}}(undef, length(data_est_other_cases) + 12),
  Vector{Union{Missing, Int64}}(undef, length(data_est_omicron_cases) + 12),
  Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations) + 12),
  vcat(data_tests, repeat(data_tests[(end - 3):end], 3))),
  posterior_samples)

CSV.write(string(results_dir, "posterior_predictive_samples_", county_id, ".csv"), DataFrame(posterior_predictive))


# MAP and MLE

# Random.seed!(1)
# map_estimate = optimize(my_model, MAP())

# map_estimate_0 = optimize(my_model, MAP(), vcat(repeat([0.0], 17), repeat([1], 2)))


# map_estimate2 = Optim.optimize(my_model, MAP(), init_vals = vcat(repeat([0.0], 17), repeat([1], 2)))


# Optim.optimize(my_model, MAP(), init_vals::AbstractArray, options::Optim.Options=Optim.Options(); kwargs...)
# function Optim.optimize(
#   model::Model, 
#   ::MAP, 
#   init_vals::AbstractArray, 
#   optimizer::Optim.AbstractOptimizer, 
#   options::Optim.Options=Optim.Options(); 
#   kwargs...
# )
# # map_estimate.f
# map_chains = Chains(transpose(reshape(repeat(map_estimate.values.array, 2000), :, 2000)),
# names(map_estimate.values,1))

# Random.seed!(1)
# map_predictive = predict(bayes_seihr(
#   Vector{Union{Missing, Int64}}(undef, length(data_est_other_cases)),
#   Vector{Union{Missing, Int64}}(undef, length(data_est_omicron_cases)),
#   Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations)),
#   vcat(data_tests)),
#   map_chains)

# CSV.write(string(results_dir, "map_predictive.csv"), DataFrame(map_predictive))

# # Random.seed!(1)
# # mle_estimate = optimize(my_model, MLE())

# # mle_chains = Chains(transpose(reshape(repeat(mle_estimate.values.array, 2000), :, 2000)),
# # names(mle_estimate.values,1))

# # Random.seed!(2)
# # mle_predictive = predict(bayes_seihr(
# #   Vector{Union{Missing, Int64}}(undef, length(data_est_other_cases) + 12),
# #   Vector{Union{Missing, Int64}}(undef, length(data_est_omicron_cases) + 12),
# #   Vector{Union{Missing, Int64}}(undef, length(data_hospitalizations) + 12),
# #   vcat(data_tests, repeat(data_tests[(end - 3):end], 3))),
# #   mle_chains)

# # CSV.write(string(results_dir, "mle_predictive.csv"), DataFrame(map_predictive))

# Array(map_predictive)[1,:]