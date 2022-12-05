prob = ODEProblem{true}(seihricud_ode_log!,
    zeros(8),
    (0.0, obstimes[end]),
    ones(9))

@model function bayes_seihricud(prob, data_new_cases, data_new_deaths, data_hospitalizations, data_icu, obstimes, param_change_times, extra_ode_precision)
    max_neg_bin_μ = popsize
    min_neg_bin_μ = 1e-4
    max_neg_bin_ϕ = 1e10
    min_neg_bin_ϕ = 1e-10

    l_obstimes = length(obstimes)
    l_param_change_times = length(param_change_times)

    # Priors
    R₀_params_non_centered ~ MvNormal(Zeros(l_param_change_times + 2), I) # +2 for sigma and init

    σ_R₀_non_centered = R₀_params_non_centered[1]
    R₀_init_non_centered = R₀_params_non_centered[2]
    log_R₀_steps_non_centered = R₀_params_non_centered[3:end]

    dur_latent_non_centered ~ Normal()
    dur_infectious_non_centered ~ Normal()
    dur_hospitalized_non_centered ~ Normal()
    dur_waning_non_centered ~ Normal()
    dur_icu_non_centered ~ Normal()

    IHR_non_centered ~ Normal()
    HICUR_non_centered ~ Normal()
    ICUDR_non_centered ~ Normal()

    case_detection_rate_non_centered ~ Normal()
    death_detection_rate_non_centered ~ Normal()

    E_init_prop_non_centered ~ Normal()
    I_init_prop_non_centered ~ Normal()
    R_init_prop_non_centered ~ Normal()

    ϕ_cases_non_centered ~ Normal()
    ϕ_hosp_non_centered ~ Normal()
    ϕ_icu_non_centered ~ Normal()
    ϕ_deaths_non_centered ~ Normal()

    # Centering tranformations
    R₀_init = exp(R₀_init_non_centered * R₀_init_non_centered_sd + R₀_init_non_centered_mean)
    σ_R₀ = exp(σ_R₀_non_centered * σ_R₀_non_centered_sd + σ_R₀_non_centered_mean)

    dur_latent = exp(dur_latent_non_centered * dur_latent_non_centered_sd + dur_latent_non_centered_mean)
    dur_infectious = exp(dur_infectious_non_centered * dur_infectious_non_centered_sd + dur_infectious_non_centered_mean)
    dur_hospitalized = exp(dur_hospitalized_non_centered * dur_hospitalized_non_centered_sd + dur_hospitalized_non_centered_mean)
    dur_waning = exp(dur_waning_non_centered * dur_waning_non_centered_sd + dur_waning_non_centered_mean)
    dur_icu = exp(dur_icu_non_centered * dur_icu_non_centered_sd + dur_icu_non_centered_mean)

    IHR = logistic(IHR_non_centered * IHR_non_centered_sd + IHR_non_centered_mean)
    HICUR = logistic(HICUR_non_centered * HICUR_non_centered_sd + HICUR_non_centered_mean)
    ICUDR = logistic(ICUDR_non_centered * ICUDR_non_centered_sd + ICUDR_non_centered_mean)

    case_detection_rate = logistic(case_detection_rate_non_centered * case_detection_rate_non_centered_sd + case_detection_rate_non_centered_mean)
    death_detection_rate = logistic(death_detection_rate_non_centered * death_detection_rate_non_centered_sd + death_detection_rate_non_centered_mean)

    E_init_prop = logistic(E_init_prop_non_centered * E_init_prop_non_centered_sd + E_init_prop_non_centered_mean)
    I_init_prop = logistic(I_init_prop_non_centered * I_init_prop_non_centered_sd + I_init_prop_non_centered_mean)
    R_init_prop = logistic(R_init_prop_non_centered * R_init_prop_non_centered_sd + R_init_prop_non_centered_mean)

    ϕ_cases = clamp.(exp(ϕ_cases_non_centered * ϕ_cases_non_centered_sd + ϕ_cases_non_centered_mean), min_neg_bin_ϕ, max_neg_bin_ϕ)
    ϕ_hosp = clamp.(exp(ϕ_hosp_non_centered * ϕ_hosp_non_centered_sd + ϕ_hosp_non_centered_mean), min_neg_bin_ϕ, max_neg_bin_ϕ)
    ϕ_icu = clamp.(exp(ϕ_icu_non_centered * ϕ_icu_non_centered_sd + ϕ_icu_non_centered_mean), min_neg_bin_ϕ, max_neg_bin_ϕ)
    ϕ_deaths = clamp.(exp(ϕ_deaths_non_centered * ϕ_deaths_non_centered_sd + ϕ_deaths_non_centered_mean), min_neg_bin_ϕ, max_neg_bin_ϕ)

    # Natural scale transformation
    γ = 1 / dur_latent
    ν = 1 / dur_infectious
    η = 1 / dur_hospitalized
    κ = 1 / dur_waning
    ω = 1 / dur_icu

    β_t_no_init = exp.(log(R₀_init) .+ cumsum(log_R₀_steps_non_centered) * σ_R₀) * ν
    β_init = R₀_init * ν

    # Initial Condititons
    E_init = E_init_prop * popsize
    I_init = I_init_prop * popsize
    R_init = R_init_prop * popsize
    # H_init is loaded as a constant
    # ICU_init is loaded as a constant
    # D_init is loaded as a constant
    # C_init is loaded as a constant
    # During optimization, E_init, I_init, R_init can sometimes get huge, such that S_init is negative
    # Should reconsider constructing this as a stick-breaking or dirichilet or some other strictly bounded prior
    S_init = max(1, remaining_population_init - (E_init + I_init + R_init))

    # ODE Setup
    function param_affect_β!(integrator)
        ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
        integrator.p[1] = β_t_no_init[ind_t] # Replace β with a new value from β_t_no_init
    end

    param_callback = PresetTimeCallback(param_change_times, param_affect_β!, save_positions=(false, false))

    abstol = extra_ode_precision ? 1e-11 : 1e-9
    reltol = extra_ode_precision ? 1e-8 : 1e-6

    u0_reg_scale = [S_init, E_init, I_init, H_init, R_init, C_init, ICU_init, D_init]
    u0_log_scale = log.(u0_reg_scale)
    p0 = [β_init, γ, ν, η, IHR, κ, HICUR, ω, ICUDR]

    sol = solve(prob, Tsit5(); callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=u0_log_scale, p=p0, tspan=(0.0, obstimes[end]))

    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end

    # Likelihood calculations
    sol_reg_scale_array = exp.(Array(sol))

    sol_hospitalizations = sol_reg_scale_array[4, 2:end]
    sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]
    sol_icu = sol_reg_scale_array[7, 2:end]
    sol_new_deaths = sol_reg_scale_array[8, 2:end] - sol_reg_scale_array[8, 1:(end-1)]

    hospitalizations_mean = clamp.(sol_hospitalizations, min_neg_bin_μ, max_neg_bin_μ)
    new_cases_mean = clamp.(sol_new_cases .* case_detection_rate, min_neg_bin_μ, max_neg_bin_μ)
    icu_mean = clamp.(sol_icu, min_neg_bin_μ, max_neg_bin_μ)
    new_deaths_mean = clamp.(sol_new_deaths .* death_detection_rate, min_neg_bin_μ, max_neg_bin_μ)

    for i in 1:l_obstimes
        data_hospitalizations[i] ~ NegativeBinomial2(hospitalizations_mean[i], ϕ_hosp)
        data_new_cases[i] ~ NegativeBinomial2(new_cases_mean[i], ϕ_cases)
        data_icu[i] ~ NegativeBinomial2(icu_mean[i], ϕ_icu)
        data_new_deaths[i] ~ NegativeBinomial2(new_deaths_mean[i], ϕ_deaths)
    end

    # Generated quantities
    S_compartment = sol_reg_scale_array[1, :]
    E_compartment = sol_reg_scale_array[2, :]
    I_compartment = sol_reg_scale_array[3, :]
    H_compartment = sol_reg_scale_array[4, :]
    R_compartment = sol_reg_scale_array[5, :]
    C_compartment = sol_reg_scale_array[6, :]
    ICU_compartment = sol_reg_scale_array[7, :]
    D_compartment = sol_reg_scale_array[8, :]

    dur_latent_days = dur_latent * time_interval_in_days
    dur_infectious_days = dur_infectious * time_interval_in_days
    dur_hospitalized_days = dur_hospitalized * time_interval_in_days
    dur_waning_days = dur_waning * time_interval_in_days
    dur_icu_days = dur_icu * time_interval_in_days

    β_t = vcat(β_init, β_t_no_init)
    R₀_t = β_t / ν
    Rₜ_t = R₀_t .* S_compartment[1:(end-1)] / popsize

    return (
        σ_R₀=σ_R₀,
        IHR=IHR,
        HICUR=HICUR,
        ICUDR=ICUDR,
        case_detection_rate=case_detection_rate,
        death_detection_rate=death_detection_rate,
        ϕ_hosp=ϕ_hosp,
        ϕ_icu=ϕ_icu,
        ϕ_cases=ϕ_cases,
        ϕ_deaths=ϕ_deaths,
        β_t=β_t,
        R₀_t=R₀_t,
        Rₜ_t=Rₜ_t,
        dur_latent_days=dur_latent_days,
        dur_infectious_days=dur_infectious_days,
        dur_hospitalized_days=dur_hospitalized_days,
        dur_waning_days=dur_waning_days,
        dur_icu_days=dur_icu_days,
        S=S_compartment,
        E=E_compartment,
        I=I_compartment,
        H=H_compartment,
        R=R_compartment,
        C=C_compartment,
        ICU=ICU_compartment,
        D=D_compartment,
        hospitalizations_mean=hospitalizations_mean,
        new_cases_mean=new_cases_mean,
        icu_mean=icu_mean,
        new_deaths_mean=new_deaths_mean
    )
end