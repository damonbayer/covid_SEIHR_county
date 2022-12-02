function seihricud_ode_log!(du, u, p, t)
    (S, E, I, H, R, C, ICU, D) = exp.(u)
    (β, γ, ν, η, IHR, κ, HICUR, ω, ICUDR) = p
    N = S + E + I + H + R + ICU + D

    # S -> E
    infection = β * I * S / N
    # E -> I
    progression = γ * E
    # I -> H
    hospitalization = ν * IHR * I
    # I -> R
    non_hospitalized_recovery = ν * (1 - IHR) * I
    # H -> ICU
    icu_admission = η * HICUR * H
    # H -> R
    hospitalized_recovery = η * (1 - HICUR) * H
    # ICU -> D
    icu_death = ω * ICUDR * ICU
    # ICU -> R
    icu_recovery = ω * (1 - ICUDR) * ICU
    # R -> S
    waned_immunity = κ * R

    @inbounds begin
        du[1] = (waned_immunity - infection) / S # S
        du[2] = (infection - progression) / E # E
        du[3] = (progression - (hospitalization + non_hospitalized_recovery)) / I # I
        du[4] = (hospitalization - (hospitalized_recovery + icu_admission)) / H # H
        du[5] = (non_hospitalized_recovery +  hospitalized_recovery + icu_recovery - waned_immunity) / R # R
        du[6] = progression / C # C
        du[7] = (icu_admission - (icu_death + icu_recovery)) / ICU # ICU
        du[8] = icu_death / D # D
    end
    nothing
end