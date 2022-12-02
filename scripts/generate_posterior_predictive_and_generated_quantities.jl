using DrWatson
using Revise
using CSV
using DataFrames
using Turing
using LogExpFunctions
using DifferentialEquations
using LinearAlgebra
using FillArrays
using Random
using covid_SEIHR_county

n_forecast_times = 12

mkpath(resultsdir("prior_generated_quantities"))
mkpath(resultsdir("posterior_generated_quantities"))
mkpath(resultsdir("prior_predictive"))
mkpath(resultsdir("posterior_predictive"))

county_id = length(ARGS) == 0 ? 28 : parse(Int64, ARGS[1])

savename_dict = Dict(:county_id => county_id)

prior_samples = load(resultsdir("prior_samples", savename("prior_samples", savename_dict, "jld2")))["prior_samples"]
posterior_samples = load(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")))["posterior_samples"]

## Load Data
include(projectdir("src/load_process_data.jl"))

## Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define ODE
include(projectdir("src/seihricud_ode_log.jl"))

## Load Model
include(projectdir("src/bayes_seihricud.jl"))

my_model = bayes_seihricud(prob, data_est_new_cases, data_est_new_deaths, data_hospitalizations, data_icu, obstimes, param_change_times, true)
my_model_forecast = bayes_seihricud(prob, data_est_new_cases_forecast, data_est_new_deaths_forecast, data_hospitalizations_forecast, data_icu_forecast, obstimes_forecast, param_change_times_forecast, true)
my_model_forecast_missing = bayes_seihricud(prob, missing_est_new_cases_forecast, missing_est_new_deaths_forecast, missing_hospitalizations_forecast, missing_icu_forecast, obstimes_forecast, param_change_times_forecast, true)

## Augment samples for forecasting
augmented_prior_samples = augment_chains_with_forecast_samples(prior_samples, my_model, my_model_forecast, "zeros")
augmented_posterior_samples = augment_chains_with_forecast_samples(posterior_samples, my_model, my_model_forecast, "zeros")

## Predictive
Random.seed!(1)
prior_predictive = predict(my_model_forecast_missing, augmented_prior_samples)
CSV.write(resultsdir("prior_predictive", savename("prior_predictive", savename_dict, "csv")), DataFrame(prior_predictive))

Random.seed!(1)
posterior_predictive = predict(my_model_forecast_missing, augmented_posterior_samples)
CSV.write(resultsdir("posterior_predictive", savename("posterior_predictive", savename_dict, "csv")), DataFrame(posterior_predictive))

## Generated Quantities
prior_generated_quantities = Chains(generated_quantities(my_model_forecast, augmented_prior_samples))
CSV.write(resultsdir("prior_generated_quantities", savename("prior_generated_quantities", savename_dict, "csv")), DataFrame(prior_generated_quantities))

posterior_generated_quantities = Chains(generated_quantities(my_model_forecast, augmented_posterior_samples))
CSV.write(resultsdir("posterior_generated_quantities", savename("posterior_generated_quantities", savename_dict, "csv")), DataFrame(posterior_generated_quantities))