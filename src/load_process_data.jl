if !@isdefined(n_forecast_times)
    @warn("n_forecast_times is not defined. Assigning n_forecast_times = 0")
    n_forecast_times = 0
end

county = subset(CSV.read("data/county_id_key.csv", DataFrame), :id => ByRow(x -> x == county_id))[1, :county]

dat = subset(CSV.read("data/cases_hospitalizations_by_county.csv", DataFrame), :county => ByRow(x -> x == county))

initialization_values = subset(CSV.read("data/initialization_values.csv", DataFrame), :county => ByRow(x -> x == county))
popsize = subset(CSV.read("data/county_pop.csv", DataFrame), :County => ByRow(x -> x == county))[!, :Population][1]

news_cases_omicron_initial = initialization_values[1, :est_omicorn_cases]
new_cases_other_initial = initialization_values[1, :est_other_cases]
hospitalizations_initial = initialization_values[1, :hospitalizations]
icu_initial = initialization_values[1, :icu]

obstimes = float(dat[:, :time])
obstimes_forecast = vcat(obstimes, obstimes[end] .+ float.(1:n_forecast_times))

param_change_times = obstimes[1:(end - 1)]
param_change_times_forecast = obstimes_forecast[1:(end - 1)]

data_est_other_cases = dat[:, :est_other_cases]
data_est_other_cases_forecast = vcat(data_est_other_cases, repeat([data_est_other_cases[end]], n_forecast_times))
missing_est_other_cases_forecast = repeat([missing], length(data_est_other_cases) + n_forecast_times)

data_est_omicron_cases = dat[:, :est_omicron_cases]
data_est_omicron_cases_forecast = vcat(data_est_omicron_cases, repeat([data_est_omicron_cases[end]], n_forecast_times))
missing_est_omicron_cases_forecast = repeat([missing], length(data_est_omicron_cases) + n_forecast_times)

data_tests = dat[:, :tests]
data_tests_forecast = vcat(data_tests, repeat([data_tests[end]], n_forecast_times))

data_est_other_tests = dat[:, :est_other_tests]
data_est_other_tests_forecast = vcat(data_est_other_tests, repeat([data_est_other_tests[end]], n_forecast_times))
end_other_cases_time = findfirst(iszero, data_est_other_cases)

if isnothing(end_other_cases_time)
  end_other_cases_time = length(data_est_other_cases) + 1
end

data_est_omicron_tests = dat[:, :est_omicron_tests]
data_est_omicron_tests_forecast = vcat(data_est_omicron_tests, repeat([data_est_omicron_tests[end]], n_forecast_times))

data_hospitalizations = dat[:, :hospitalizations]
data_hospitalizations_forecast = vcat(data_hospitalizations, repeat([data_hospitalizations[end]], n_forecast_times))
missing_hospitalizations_forecast = repeat([missing], length(data_hospitalizations) + n_forecast_times)

data_icu = dat[:, :icu]
data_icu_forecast = vcat(data_icu, repeat([data_icu[end]], n_forecast_times))
missing_icu_forecast = repeat([missing], length(data_icu) + n_forecast_times)

data_est_death = dat[:, :est_deaths]
data_est_death_forecast = vcat(data_est_death, repeat([data_est_death[end]], n_forecast_times))
missing_est_death_forecast = repeat([missing], length(data_est_death) + n_forecast_times)

