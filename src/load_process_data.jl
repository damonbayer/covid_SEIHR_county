if !@isdefined(n_forecast_times)
    @warn("n_forecast_times is not defined. Assigning n_forecast_times = 0")
    n_forecast_times = 0
end

county = subset(CSV.read("data/county_id_pop.csv", DataFrame), :id => ByRow(x -> x == county_id))[1, :county]

dat = subset(CSV.read("data/cases_hospitalizations_by_county.csv", DataFrame), :county => ByRow(x -> x == county))
initialization_values = subset(CSV.read("data/initialization_values.csv", DataFrame), :county => ByRow(x -> x == county))

const popsize = float(initialization_values[1, :population])
const init_H = float(initialization_values[1, :H])
const init_ICU = float(initialization_values[1, :ICU])
const init_D = float(initialization_values[1, :D])
const init_remaining_population = float(initialization_values[1, :remaining_population])

obstimes = float(dat[:, :time])
obstimes_forecast = vcat(obstimes, obstimes[end] .+ float.(1:n_forecast_times))

param_change_times = obstimes[1:(end-1)]
param_change_times_forecast = obstimes_forecast[1:(end-1)]

data_est_cases = dat[:, :est_cases]
data_est_cases_forecast = vcat(data_est_cases, repeat([data_est_cases[end]], n_forecast_times))
missing_est_cases_forecast = repeat([missing], length(data_est_cases) + n_forecast_times)

data_hospitalizations = dat[:, :hospitalized]
data_hospitalizations_forecast = vcat(data_hospitalizations, repeat([data_hospitalizations[end]], n_forecast_times))
missing_hospitalizations_forecast = repeat([missing], length(data_hospitalizations) + n_forecast_times)

data_icu = dat[:, :icu]
data_icu_forecast = vcat(data_icu, repeat([data_icu[end]], n_forecast_times))
missing_icu_forecast = repeat([missing], length(data_icu) + n_forecast_times)

data_est_death = dat[:, :est_deaths]
data_est_death_forecast = vcat(data_est_death, repeat([data_est_death[end]], n_forecast_times))
missing_est_death_forecast = repeat([missing], length(data_est_death) + n_forecast_times)