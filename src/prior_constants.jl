const other_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.1)
const omicron_case_detection_rate_samples = logistic.(randn(2000) * 0.2 .- 1.6)
const other_rho_samples = other_case_detection_rate_samples / median(data_est_other_tests)
const omicron_rho_samples = omicron_case_detection_rate_samples / median(data_est_omicron_tests)

const other_rho_mean = mean(logit.(other_rho_samples))
const other_rho_sd = std(logit.(other_rho_samples))

const omicron_rho_mean = mean(logit.(omicron_rho_samples))
const omicron_rho_sd = std(logit.(omicron_rho_samples))