const R₀_init_non_centered_mean = 0.55
const R₀_init_non_centered_sd = 0.2

const σ_R0_non_centered_mean = -3
const σ_R0_non_centered_sd = 0.2

const case_detection_rate_non_centered_mean = -1.78
const case_detection_rate_non_centered_sd = 0.1

const dur_latent_non_centered_mean = -1.5
const dur_latent_non_centered_sd = 0.25

const dur_infectious_non_centered_mean = log(5 / 7)
const dur_infectious_non_centered_sd = 0.25

const IHR_non_centered_mean = -4.3
const IHR_non_centered_sd = 0.25

const dur_hospitalized_non_centered_mean = -1.54
const dur_hospitalized_non_centered_sd = 0.1

const dur_waning_non_centered_mean = log(12)
const dur_waning_non_centered_sd = 0.1

const dur_icu_non_centered_mean = log(0.31)
const dur_icu_non_centered_sd = 0.1

const HICUR_non_centered_mean = -1.69
const HICUR_non_centered_sd = 0.2

const ICUDR_non_centered_mean = -1.59
const ICUDR_non_centered_sd = 0.2

const death_detection_rate_non_centered_mean = 2.3
const death_detection_rate_non_centered_sd = 0.2

overdisp_priors = CSV.read(datadir(string("overdisp_priors/overdisp_priors_countyid=", county_id, ".csv")), DataFrame)
const ϕ_hosp_non_centered_sd = overdisp_priors[overdisp_priors.datastream .== "hospitalized", :sd][1]
const ϕ_hosp_non_centered_mean = overdisp_priors[overdisp_priors.datastream .== "hospitalized", :mean][1]
const ϕ_icu_non_centered_sd = overdisp_priors[overdisp_priors.datastream .== "icu", :sd][1]
const ϕ_icu_non_centered_mean = overdisp_priors[overdisp_priors.datastream .== "icu", :mean][1]
