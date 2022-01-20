# covid_SEIHR_county

This repository contains a simple model used for forecasting hospitalizations in California counties in January 2022.

# Data

* [Hospitalization](https://data.ca.gov/dataset/covid-19-hospital-data1) and [case](https://data.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state1) data come from the [California Open Data Portal](https://data.ca.gov).
* Delays for case reporting are estimated using line lists provided by the [Orange County Health Care Agency](https://www.ochealthinfo.com).
* [Prevalence of variants in California](https://github.com/blab/rt-from-frequency-dynamics/blob/master/data/omicron-us/omicron-us_location-variant-sequence-counts.tsv) is estimated from from [GISAID](https://www.gisaid.org) data via the [Bedford Lab](https://bedford.io)

We use `est_omicron_cases = reported_cases / proportion_cases_expected * proportion_omicron` and `est_other_cases = reported_cases / proportion_cases_expected * (1 - proportion_omicron)` as data in the model.

# Model

The model is a S-E-I-R type model startified by omicron vs non-omicrons  with hospitalization compartments.
The model is fit with [fit_model.jl](fit_model.jl) using [Turing.jl](https://turing.ml/stable/).

![model_diagram](model_diagram.png)

Priors and posteriors for the model parameters are presented in [prior_post_plots.pdf](prior_post_plots.pdf).

# Results

Results for each county are presented in [post_pred_plots.pdf](post_pred_plots.pdf).

