# distribution notes form ICU and Death model


# previously used hospitalization duration priors -------------------------


  dur_hosp_other <- exp(rnorm(10000) * 0.1 - 0.36)
  quantile(dur_hosp_other, c(0.025, 0.5, 0.975)) * 7

  dur_hosp_omi <- exp(rnorm(10000) * 0.1 - 1.54)
  quantile(dur_hosp_omi, c(0.025, 0.5, 0.975)) * 7

  # the mean of this distribution comes from LEMMA
  # https://localepi.github.io/LEMMA/ last accessed on 9/12/2022
  # SD comes from previously used SD for dur_hosp_other
  dur_icu_other <- exp(rnorm(10000) * 0.1 + log(7/7))
  quantile(dur_icu_other, c(0.025, 0.5, 0.975)) * 7

  # the ratio of the mean of the duration hospitalization for Omicron vs non-Omicron
  # was previously set to be exp(-1.54)/exp(0.36) = 0.31 (approx)
  # I preserved this ratio for the duration icu for omicron
  dur_icu_omi <- exp(rnorm(10000) * 0.1 + log(0.31))
  quantile(dur_icu_omi, c(0.025, 0.5, 0.975)) * 7



# percent of hospitalized who go to ICU -----------------------------------

  # the mean of this distribution is from LEMMA
  # https://localepi.github.io/LEMMA/ last accessed on 9/12/2022
  # looks like I stole the SD from previously used SD for infection to Hospitalized distribution
  HICUR_other <- plogis(rnorm(10000) * 0.2 + qlogis(0.21))
  quantile(HICUR_other, c(0.025, 0.5, 0.975))

  #cdc information come from: https://www.cdc.gov/mmwr/volumes/71/wr/pdfs/mm7104e4-H.pdf
  # cdc says rates are down about 26%, I scaled the mean of the LEMMA distribution by 0.74
  HICUR_omicron <- plogis(rnorm(10000) * 0.2 + qlogis(0.21 * 0.74))
  quantile(HICUR_omicron, c(0.025, 0.5, 0.975))

  # ICU to deaths comes from LEMMA
  ICUDR_other <- plogis(rnorm(10000) * 0.2 + qlogis(0.25))
  quantile(ICUDR_other, c(0.025, 0.5, 0.975))

  # cdc says death in hospital is down 42%, I scaled the mean of the LEMMA distr. by 0.68
  ICUDR_omicron <- plogis(rnorm(10000) * 0.2 + qlogis(0.25 * 0.68))
  quantile(ICUDR_omicron, c(0.025, 0.5, 0.975))


# previous hospitalization fraction priors ---------------------------------
  IHR_other <- plogis(rnorm(10000) * 0.2 - 3.2)
  quantile(IHR_other)



  IHR_omicron <- plogis(rnorm(10000) * 0.25 - 4.3)
  quantile(IHR_omicron)

# looking at overall infection to fatality ratio --------------------------
  IFR_other <- IHR_other * HICUR_other * ICUDR_other
  IFR_omicron <- IHR_omicron * HICUR_omicron * ICUDR_omicron

  quantile(IFR_other)
  quantile(IFR_omicron)

  median(IHR_other) * median(HICUR_other) * median(ICUDR_other)
  median(IHR_omicron) * median(HICUR_omicron) * median(ICUDR_omicron)



# death detection rate ----------------------------------------------------
  # looks like I took it from Damon's manuscript
  # although his currently lists Logit-Normal(2.3, 0.04) (as of 9/12/2022)
  # However, our quantiles match perfectly, so I still think this is
  # where the parameters came from
  death_detection_rate <- plogis(rnorm(10000) * 0.2 + 2.313635)
  quantile(death_detection_rate, c(0.025, 0.5, 0.975))


