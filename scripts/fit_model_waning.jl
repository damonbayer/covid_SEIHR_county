using DrWatson
using Revise
using JLD2
using FileIO
using CSV
using DataFrames
using Turing
using LinearAlgebra
using FillArrays
using DifferentialEquations
using LogExpFunctions
using Random
using ForwardDiff
using Optim
using Random
using LineSearches
using ForwardDiff
using covid_SEIHR_county

county_id = length(ARGS) == 0 ? 1 : parse(Int64, ARGS[1])

priors_only = county_id == 0

if priors_only
  county_id = 1
end

mkpath(resultsdir())
mkpath(resultsdir("posterior_samples"))

savename_dict = Dict(:county_id => county_id)

## Control Parameters
n_samples = 2_000
n_chains = 4
time_interval_in_days = 7

## Load Data
include(projectdir("src/load_process_data.jl"))

## Load overdisp priors_only
overdisp_priors = CSV.read(datadir(string("overdisp_priors/overdisp_priors_countyid=", county_id, ".csv")), DataFrame)
const ϕ_hosp_sd = overdisp_priors[overdisp_priors.datastream .== "hosp", :sd][1]
const ϕ_hosp_mean = overdisp_priors[overdisp_priors.datastream .== "hosp", :mean][1]
const ϕ_icu_sd = overdisp_priors[overdisp_priors.datastream .== "icu", :sd][1]
const ϕ_icu_mean = overdisp_priors[overdisp_priors.datastream .== "icu", :mean][1]
## Define Priors
include(projectdir("src/prior_constants.jl"))

## Define ODE
include(projectdir("src/seir_ode_log.jl"))

## Load Model
include(projectdir("src/bayes_seihr_waning.jl"))

my_model = bayes_seihr(
  prob,
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

# Sample Prior
if priors_only
  Random.seed!(county_id)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
  wsave(resultsdir("prior_samples.jld2"), @dict prior_samples)
  exit()
end

MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]

alg = Gibbs(NUTS(-1, 0.8,
    :prop_omicron_only_init_non_centered,
    :dur_latent_non_centered_non_omicron,
    :dur_infectious_non_centered_non_omicron,
    :IHR_non_centered_non_omicron,
    :dur_hospitalized_non_centered_non_omicron,
    :dur_icu_non_centered_non_omicron,
    :HICUR_non_centered_non_omicron,
    :ICUDR_non_centered_non_omicron,
    :E_init_non_centered_non_omicron,
    :I_init_non_centered_non_omicron,
    :case_detection_rate_non_centered_other,
    :dur_latent_non_centered_omicron,
    :dur_infectious_non_centered_omicron,
    :IHR_non_centered_omicron,
    :dur_hospitalized_non_centered_omicron,
    :dur_icu_non_centered_omicron,
    :HICUR_non_centered_omicron,
    :ICUDR_non_centered_omicron,
    :dur_waning_non_centered_omicron,
    :E_init_non_centered_omicron,
    :I_init_non_centered_omicron,
    :case_detection_rate_non_centered_omicron,
    :death_detection_rate_non_centered,
    :ϕ_cases_non_centered,
    :ϕ_hospitalizations_non_centered,
    :ϕ_death_non_centered,
    :ϕ_icu_non_centered),
  ESS(:R0_params_non_centered))

# Sample Posterior
Random.seed!(county_id)
posterior_samples = sample(my_model, alg, MCMCThreads(), n_samples, n_chains, init_params = repeat([MAP_init], n_chains) .* collect(range(0.92, stop = 0.98, length = n_chains)))

wsave(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")), @dict posterior_samples)
