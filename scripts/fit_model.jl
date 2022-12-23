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

mkpath(resultsdir())
mkpath(resultsdir("prior_samples"))
mkpath(resultsdir("posterior_samples"))

county_id = length(ARGS) == 0 ? 1 : parse(Int64, ARGS[1])
savename_dict = Dict(:county_id => county_id)

## Control Parameters
n_samples = 250
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
MAP_values = optimize_many_MAP(my_model, 10, 1, true)[1]
Random.seed!(1)
MAP_noise = [randn(length(MAP_values)) for _ in 1:n_chains]
init_params = repeat([MAP_values], n_chains) * 0.95 + MAP_noise * 0.05

alg = NUTS(-1, 0.8)

Random.seed!(1)
posterior_samples = sample(my_model,
    alg,
    MCMCThreads(),
    n_samples,
    n_chains,
    init_params=init_params)

wsave(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")), @dict posterior_samples)
