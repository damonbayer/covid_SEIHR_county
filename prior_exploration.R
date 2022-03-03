expit <- logistic <- plogis
logit <- qlogis
target_qunatiles <- c(0.05, 0.5, 0.95)
z <- rnorm(10000)

# R₀ₙ
exp(qnorm(p = target_qunatiles, mean = 0.25, sd = 0.2))
# R₀ₒ
exp(qnorm(p = target_qunatiles, mean = 1.25, sd = 0.2))

exp(qnorm(p = target_qunatiles, mean = 1.5, sd = 0.2)) * 0.3

# dur_latentₙ_days
exp(qnorm(p = target_qunatiles, mean = -1.27, sd = 0.25)) * 7

# dur_latentₒ_days
exp(qnorm(p = target_qunatiles, mean = -1.5, sd = 0.25)) * 7


# https://dash.harvard.edu/handle/1/37370587
# dur_infectiousₙ_days
exp(qnorm(p = target_qunatiles, mean = -0.568, sd = 0.25)) * 7
# dur_infectiousₒ
exp(qnorm(p = target_qunatiles, mean = -0.85, sd = 0.25)) * 7

# https://twitter.com/famulare_mike/status/1482097499045195779/photo/1
# IHRₙ
scales::percent(expit(qnorm(p = target_qunatiles, mean = -4.2, sd = 0.2)))
scales::percent(expit(qnorm(p = target_qunatiles, mean = -3.3, sd = 0.25)))
# IHRₒ
scales::percent(expit(qnorm(p = target_qunatiles, mean = -5.3, sd = 0.25)))
scales::percent(expit(qnorm(p = target_qunatiles, mean = -4.3, sd = 0.25)))


# https://www.medrxiv.org/content/10.1101/2022.01.11.22269045v1

# dur_hospitalizedₙ
exp(qnorm(p = target_qunatiles, mean = -0.36, sd = 0.1)) * 7
# dur_hospitalizedₒ
exp(qnorm(p = target_qunatiles, mean = -1.54, sd = 0.1)) * 7



1/expit(qnorm(p = target_qunatiles, mean = -1.1, sd = 0.2))

logistic.(randn(2000) * 0.2 .- 1.1)



# prop_omicron_only
expit(qnorm(p = target_qunatiles, mean = 1.2, sd = 0.15))
