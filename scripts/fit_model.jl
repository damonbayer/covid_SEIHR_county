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
using Optim
using LineSearches
using covid_SEIHR_county

mkpath(resultsdir("prior_samples"))
mkpath(resultsdir("posterior_samples"))

county_id = length(ARGS) == 0 ? 1 : parse(Int64, ARGS[1])

mkpath(resultsdir())
mkpath(resultsdir("posterior_samples"))

savename_dict = Dict(:county_id => county_id)

## Control Parameters
n_samples = 2_000
n_chains = 4

## Load Data
include(projectdir("src/load_process_data.jl"))

## Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define ODE
include(projectdir("src/seihricud_ode_log.jl"))

## Load Model
include(projectdir("src/bayes_seihricud.jl"))

my_model = bayes_seihricud(prob, data_est_new_cases, data_est_new_deaths, data_hospitalizations, data_icu, obstimes, param_change_times, false)

## Sample Prior
Random.seed!(1)
prior_samples = sample(my_model, Prior(), 2_000)
wsave(resultsdir("prior_samples", savename("prior_samples", savename_dict, "jld2")), @dict prior_samples)

## Fit Posterior
MAP_values = optimize_many_MAP(my_model, 20, 1, true)[1]
Random.seed!(1)
MAP_noise = [vcat(randn(length(MAP_values) - 2), rand(Exponential(), 2)) for _ in 1:n_chains]
init_params = repeat([MAP_values], n_chains) * 0.95 + MAP_noise * 0.05

alg = Gibbs(NUTS(-1, 0.8,
        :dur_latent_non_centered,
        :dur_infectious_non_centered,
        :dur_hospitalized_non_centered,
        :dur_waning_non_centered,
        :dur_icu_non_centered,
        :IHR_non_centered,
        :HICUR_non_centered,
        :ICUDR_non_centered,
        :case_detection_rate_non_centered,
        :death_detection_rate_non_centered,
        :E_init_prop_non_centered,
        :I_init_prop_non_centered,
        :R_init_prop_non_centered,
        :ϕ_hosp_non_centered,
        :ϕ_icu_non_centered,
        :ϕ_cases_non_centered,
        :ϕ_deaths_non_centered),
    ESS(:R₀_params_non_centered))

Random.seed!(1)
posterior_samples = sample(my_model,
    alg,
    MCMCThreads(),
    n_samples,
    n_chains,
    init_params=init_params)

wsave(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")), @dict posterior_samples)
