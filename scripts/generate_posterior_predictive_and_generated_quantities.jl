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
  0
else
  parse(Int64, ARGS[1])
end

priors_only = county_id == 0

if priors_only
  county_id = 1
end

mkpath(resultsdir("generated_quantities"))
mkpath(resultsdir("posterior_predictive"))

savename_dict = Dict(:county_id => county_id)

## Control Parameters
n_forecast_times = 12

## Load Data
include(projectdir("src/load_process_data.jl"))

## Load overdisp priors_only
overdisp_priors = CSV.read(datadir(string("overdisp_priors_countyid", county_id, ".csv")), DataFrame)
ϕ_hosp_sd = overdisp_priors[overdisp_priors.labels .== "hosp", :sd] 
ϕ_hosp_mean = overdisp_priors[overdisp_priors.labels .== "hosp", :mean]
ϕ_icu_sd = overdisp_priors[overdisp_priors.labels .== "icu", :sd]
ϕ_icu_mean = overdisp_priors[overdisp_priors.labels .== "icu", :mean]

## Define Priors
include(projectdir("src/prior_constants.jl"))

## Define ODE
include(projectdir("src/seir_ode_log.jl"))

## Load Model
include(projectdir("src/bayes_seihr_waning.jl"))

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
  missing_icu_forecast, 
  missing_est_death_forecast,
  obstimes_forecast,
  param_change_times_forecast,
  true,
  end_other_cases_time)

if priors_only
    prior_samples = load(resultsdir("prior_samples.jld2"))["prior_samples"]
    
    Random.seed!(county_id)
    prior_samples_forecast_zeros = augment_chains_with_forecast_samples(Chains(prior_samples, :parameters), my_model, my_model_forecast, "zeros")
    prior_indices_to_keep = .!isnothing.(generated_quantities(my_model_forecast, prior_samples_forecast_zeros));
    
    Random.seed!(county_id)
    prior_predictive_zeros = predict(my_model_forecast_missing, prior_samples_forecast_zeros)
    CSV.write(resultsdir("prior_predictive.csv"), DataFrame(prior_predictive_zeros))
    
    Random.seed!(county_id)
    gq_zeros = get_gq_chains(my_model_forecast, prior_samples_forecast_zeros);
    CSV.write(resultsdir("prior_generated_quantities.csv"), DataFrame(gq_zeros))
    exit()
end

posterior_samples = load(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")))["posterior_samples"]

Random.seed!(county_id)
posterior_samples_forecast_zeros = augment_chains_with_forecast_samples(Chains(posterior_samples, :parameters), my_model, my_model_forecast, "zeros")
indices_to_keep = .!isnothing.(generated_quantities(my_model_forecast, posterior_samples_forecast_zeros));
posterior_samples_forecast_zeros = ChainsCustomIndex(posterior_samples_forecast_zeros, indices_to_keep);

Random.seed!(county_id)
predictive_zeros = predict(my_model_forecast_missing, posterior_samples_forecast_zeros)
CSV.write(resultsdir("posterior_predictive", savename("posterior_predictive", savename_dict, "csv")), DataFrame(predictive_zeros))

Random.seed!(county_id)
gq_zeros = get_gq_chains(my_model_forecast, posterior_samples_forecast_zeros);
CSV.write(resultsdir("generated_quantities", savename("generated_quantities", savename_dict, "csv")), DataFrame(gq_zeros))
