# =============================================================================
# Task B: Synthetic Data Generation
# =============================================================================
#
# Generate synthetic TGA (Thermogravimetric Analysis) data for all 4 experiments
# using the true parameters. Add Gaussian noise to simulate real sensor data.
# Output a DataFrame in PEtab-compatible format.
# =============================================================================

using OrdinaryDiffEq, DataFrames, Random

"""
    simulate_experiment(rn, exp_id, exp_cond, params; n_points=150)

Simulate a single TGA experiment using the Catalyst.jl reaction system.

Arguments:
- `rn`: The completed ReactionSystem
- `exp_id`: Experiment identifier string (e.g., "Exp1")
- `exp_cond`: Named tuple with (beta, PO2) for this experiment
- `params`: Dict of true parameter values
- `n_points`: Number of data points to sample

Returns: (times, temperatures, mass_values) arrays
"""
function simulate_experiment(rn, exp_id, exp_cond, params; n_points=150)
    # Compute simulation time span from temperature range
    # T = T0 + β*t  →  t_final = (T_final - T0) / β
    t_final = (T_FINAL - INITIAL_TEMP) / exp_cond.beta

    # Map species to initial conditions
    @unpack M, C, F, Temp, G1, G2 = rn
    u0 = [
        M    => INITIAL_M,
        C    => INITIAL_C,
        F    => INITIAL_F,
        Temp => INITIAL_TEMP,
        G1   => INITIAL_G1,
        G2   => INITIAL_G2,
    ]

    # Map all parameters (true values + experimental conditions)
    @unpack log10_A1, E1, n1, nu_char, log10_A2, E2, n2, m2,
            log10_A3, E3, n3, m3, PO2, beta = rn
    pmap = [
        log10_A1 => params[:log10_A1],
        E1       => params[:E1],
        n1       => params[:n1],
        nu_char  => params[:nu_char],
        log10_A2 => params[:log10_A2],
        E2       => params[:E2],
        n2       => params[:n2],
        m2       => params[:m2],
        log10_A3 => params[:log10_A3],
        E3       => params[:E3],
        n3       => params[:n3],
        m3       => params[:m3],
        PO2      => exp_cond.PO2,
        beta     => exp_cond.beta,
    ]

    # Create and solve the ODE problem
    oprob = ODEProblem(rn, u0, (0.0, t_final), pmap)
    sol = solve(oprob, Rodas5P(); abstol=1e-10, reltol=1e-8, saveat=t_final/n_points)

    # Extract results
    times = sol.t
    # Mass = M + C + F (solid mass fraction — excludes gases)
    mass_values = sol[M] .+ sol[C] .+ sol[F]
    temperatures = sol[Temp]

    return times, temperatures, mass_values
end

"""
    generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)

Generate noisy synthetic TGA data for all 4 experiments.

Arguments:
- `rn`: The completed ReactionSystem
- `noise_sigma`: Gaussian noise standard deviation (0.5% = 0.005)
- `n_points`: Number of data points per experiment
- `seed`: Random seed for reproducibility

Returns: (measurements_df, clean_data)
- `measurements_df`: PEtab-compatible DataFrame
- `clean_data`: Dict of clean simulation results for each experiment
"""
function generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)
    Random.seed!(seed)

    # Storage for clean data (for plotting reference)
    clean_data = Dict{String, NamedTuple}()

    # PEtab measurement DataFrame
    all_obs_ids    = String[]
    all_sim_ids    = String[]
    all_times      = Float64[]
    all_measurements = Float64[]

    for (exp_id, exp_cond) in sort(collect(EXPERIMENTS))
        println("  Simulating $exp_id: β = $(exp_cond.beta*60) K/min, PO₂ = $(exp_cond.PO2)")

        times, temperatures, mass_clean = simulate_experiment(
            rn, exp_id, exp_cond, TRUE_PARAMS; n_points=n_points
        )

        # Add Gaussian noise (σ = 0.5% of initial mass = 0.005)
        noise = noise_sigma .* randn(length(mass_clean))
        mass_noisy = mass_clean .+ noise
        # Clamp to physically reasonable range [0, 1]
        mass_noisy = clamp.(mass_noisy, 0.0, 1.0)

        # Store clean data for later plotting
        clean_data[exp_id] = (
            times = times,
            temperatures = temperatures,
            mass_clean = mass_clean,
            mass_noisy = mass_noisy,
        )

        # Append to PEtab measurement arrays
        for i in eachindex(times)
            push!(all_obs_ids, "mass_total")
            push!(all_sim_ids, exp_id)
            push!(all_times, times[i])
            push!(all_measurements, mass_noisy[i])
        end
    end

    # Build PEtab-compatible measurements DataFrame
    measurements_df = DataFrame(
        obs_id        = all_obs_ids,
        simulation_id = all_sim_ids,
        time          = all_times,
        measurement   = all_measurements,
    )

    return measurements_df, clean_data
end
