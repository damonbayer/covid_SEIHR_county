function seir_ode_log!(du, u, p, t)
    (S_both, S_omicron_only, E_non_omicron, E_omicron, I_non_omicron, I_omicron, H_non_omicron, H_omicron, R, C_non_omicron, C_omicron, 
    ICU_non_omicron, 
    ICU_omicron,
    D_non_omicron,
    D_omicron) = exp.(u)
    (β_init_non_omicron, β_init_omicron, γ_non_omicron, γ_omicron, ν_non_omicron, ν_omicron, η_non_omicron, η_omicron, IHR_non_omicron, IHR_omicron, κ,
    HICUR_non_omicron, HICUR_omicron, ω_non_omicron, ω_omicron, ICUDR_non_omicron, ICUDR_omicron) = p
    N = S_both + S_omicron_only + E_non_omicron + E_omicron + I_non_omicron + I_omicron + H_non_omicron + H_omicron + R + ICU_non_omicron + ICU_omicron + D_omicron + D_non_omicron

  # infection = β * I * S / N
  infection_omicron_only_to_omicron = β_init_omicron * I_omicron * S_omicron_only / N
  infection_both_to_omicron = β_init_omicron * I_omicron * S_both / N
  infection_both_to_non_omicron = β_init_non_omicron * I_non_omicron * S_both / N
  
  infection_omicron = infection_both_to_omicron + infection_omicron_only_to_omicron
  infection_non_omicron = infection_both_to_non_omicron

  # progression = γ * E
  progression_non_omicron = γ_non_omicron * E_non_omicron
  progression_omicron = γ_omicron * E_omicron

  # hospitalization = ν * IHR * I
  hospitalization_non_omicron = ν_non_omicron * IHR_non_omicron * I_non_omicron
  hospitalization_omicron = ν_omicron * IHR_omicron * I_omicron

  # non_hospitalized_recovery = ν * (1 - IHR) * I
  non_hospitalized_recovery_non_omicron = ν_non_omicron * (1 - IHR_non_omicron) * I_non_omicron
  non_hospitalized_recovery_omicron = ν_omicron * (1 - IHR_omicron) * I_omicron

  # hospitalized_recovery = η * H
  hospitalized_recovery_non_omicron = η_non_omicron * (1-HICUR_non_omicron) * H_non_omicron
  hospitalized_recovery_omicron = η_omicron * (1-HICUR_omicron) * H_omicron

  # waning immunity = \kappa * R, and I will assume it only goes back to S_omicron_only
  waning_immunity = κ * R

  # progression from H to ICU_non_omicron
  icu_non_omicron = η_non_omicron * HICUR_non_omicron * H_non_omicron
  icu_omicron = η_omicron * HICUR_omicron * H_omicron

  # progression from ICU to deaths 
  death_non_omicron = ω_non_omicron * ICUDR_non_omicron * ICU_non_omicron
  death_omicron = ω_omicron * ICUDR_omicron * ICU_omicron 

  # progression from ICU to recovery 
  icu_recovery_non_omicron = ω_non_omicron * (1-ICUDR_non_omicron) * ICU_non_omicron
  icu_recovery_omicron = ω_omicron * (1-ICUDR_omicron) * ICU_omicron



  @inbounds begin
    du[1] = -(infection_both_to_omicron + infection_both_to_non_omicron) / S_both # S_both
    du[2] = (waning_immunity - infection_omicron_only_to_omicron) / S_omicron_only # S_omicron_only
    du[3] = (infection_non_omicron - progression_non_omicron) / E_non_omicron # E_non_omicron
    du[4] = (infection_omicron - progression_omicron) / E_omicron # E_omicron
    du[5] = (progression_non_omicron - (hospitalization_non_omicron + non_hospitalized_recovery_non_omicron)) / I_non_omicron # I_non_omicron
    du[6] = (progression_omicron - (hospitalization_omicron + non_hospitalized_recovery_omicron)) / I_omicron # I_omicron
    du[7] = (hospitalization_non_omicron - hospitalized_recovery_non_omicron) / H_non_omicron # H_non_omicron
    du[8] = (hospitalization_omicron - hospitalized_recovery_omicron) / H_omicron # H_omicron
    du[9] = (non_hospitalized_recovery_non_omicron + hospitalized_recovery_non_omicron + non_hospitalized_recovery_omicron + hospitalized_recovery_omicron + icu_recovery_non_omicron + 
    icu_recovery_omicron - waning_immunity) / R  # R
    du[10] = progression_non_omicron / C_non_omicron # C_non_omicron
    du[11] = progression_omicron / C_omicron # C_omicron
    du[12] = (icu_non_omicron - icu_recovery_non_omicron - death_non_omicron)/ICU_non_omicron
    du[13] = (icu_omicron - icu_recovery_omicron - death_omicron)/ICU_omicron
    du[14] = (death_non_omicron)/D_non_omicron
    du[15] = (death_omicron)/D_omicron
  end
    nothing
end