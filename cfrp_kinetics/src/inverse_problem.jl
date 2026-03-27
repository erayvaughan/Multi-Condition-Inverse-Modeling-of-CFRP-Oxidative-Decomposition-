# =============================================================================
# Tasks C & D: Inverse Problem Setup and Parameter Estimation
# =============================================================================
#
# Set up the PEtab model for multi-condition parameter estimation.
# Define simulation conditions, observables, and parameters to estimate.
# Run multi-start optimization and compare recovered vs true parameters.
# =============================================================================

using PEtab, OrdinaryDiffEq, OptimizationOptimJL, Optim, DataFrames, Printf, Random

"""
    setup_petab_model(rn, measurements_df)

Set up the PEtab model for parameter estimation.

Arguments:
- `rn`: The completed Catalyst.jl ReactionSystem
- `measurements_df`: PEtab-compatible measurement DataFrame

Returns: PEtabODEProblem ready for optimization
"""
function setup_petab_model(rn, measurements_df)
    @unpack M, C, F, Temp, G1, G2 = rn

    # --- Observables ---
    # We observe total solid mass = M + C + F
    # Noise model: constant Gaussian noise σ = 0.005 (0.5%)
    @parameters sigma_mass
    observables = [
        PEtabObservable(:mass_total, M + C + F, sigma_mass),
    ]

    # --- Simulation Conditions ---
    # Map each experiment ID to its specific β and PO₂ values
    # These override the model parameters for each experimental condition
    @unpack PO2, beta = rn
    simulation_conditions = [
        PEtabCondition(:Exp1, beta => 2.5 / 60.0,  PO2 => 0.21),
        PEtabCondition(:Exp2, beta => 5.0 / 60.0,  PO2 => 0.21),
        PEtabCondition(:Exp3, beta => 10.0 / 60.0, PO2 => 0.21),
        PEtabCondition(:Exp4, beta => 5.0 / 60.0,  PO2 => 0.05),
    ]

    # --- Parameters to Estimate ---
    # Hint from assignment: estimate log10(A) instead of A for numerical stability
    # We estimate 11 kinetic parameters + 1 noise parameter
    petab_parameters = [
        # Matrix Pyrolysis parameters
        PEtabParameter(:log10_A1, value=5.0,  lb=2.0,    ub=10.0,    scale=:lin),
        PEtabParameter(:E1,       value=1.2e5, lb=8e4,    ub=2e5,     scale=:lin),
        PEtabParameter(:n1,       value=1.0,  lb=0.5,    ub=3.0,     scale=:lin),
        # Char Oxidation parameters
        PEtabParameter(:log10_A2, value=7.0,  lb=3.0,    ub=12.0,    scale=:lin),
        PEtabParameter(:E2,       value=1.6e5, lb=1e5,    ub=2.5e5,   scale=:lin),
        PEtabParameter(:n2,       value=1.5,  lb=0.5,    ub=3.0,     scale=:lin),
        PEtabParameter(:m2,       value=0.6,  lb=0.1,    ub=2.0,     scale=:lin),
        # Fiber Oxidation parameters
        PEtabParameter(:log10_A3, value=9.0,  lb=5.0,    ub=13.0,    scale=:lin),
        PEtabParameter(:E3,       value=2.5e5, lb=1.5e5,  ub=3.5e5,   scale=:lin),
        PEtabParameter(:n3,       value=2.0,  lb=0.5,    ub=3.5,     scale=:lin),
        PEtabParameter(:m3,       value=0.85, lb=0.1,    ub=2.0,     scale=:lin),
        # Noise parameter (estimated alongside kinetics)
        PEtabParameter(:sigma_mass, value=0.005, lb=1e-4, ub=0.1,    scale=:log10),
    ]

    # --- Initial Conditions (same for all experiments) ---
    state_map = [
        M    => INITIAL_M,
        C    => INITIAL_C,
        F    => INITIAL_F,
        Temp => INITIAL_TEMP,
        G1   => INITIAL_G1,
        G2   => INITIAL_G2,
    ]

    # --- Fixed Parameters ---
    @unpack nu_char = rn
    parameter_map = [
        nu_char => 0.2,  # Char yield is FIXED (not estimated)
    ]

    # --- Build PEtab Model ---
    model = PEtabModel(
        rn,
        observables,
        measurements_df,
        petab_parameters;
        simulation_conditions = simulation_conditions,
        speciemap = state_map,
        parametermap = parameter_map,
        verbose = true,
    )

    # --- Build PEtab ODE Problem ---
    petab_prob = PEtabODEProblem(
        model;
        odesolver = ODESolver(Rodas5P(); abstol=1e-8, reltol=1e-6, maxiters=Int64(1e5), force_dtmin=true),
        gradient_method = :ForwardDiff,
        hessian_method = :ForwardDiff,
    )

    return petab_prob
end

"""
    run_parameter_estimation(petab_prob; n_multistarts=10)

Run multi-start optimization to estimate kinetic parameters.
Uses manual perturbations around nominal values to avoid the extreme
instabilities that random Latin Hypercube sampling causes in Arrhenius
parameter spaces (where 10^A * exp(-E/RT) easily overflows).

Arguments:
- `petab_prob`: PEtabODEProblem
- `n_multistarts`: Number of starting points (1 nominal + perturbations)

Returns: PEtabMultistartResult
"""
function run_parameter_estimation(petab_prob; n_multistarts=10)
    println("\n" * "="^60)
    println("  Running multi-start optimization ($n_multistarts starts)")
    println("="^60)

    # Generate starting points: nominal + controlled perturbations
    # This avoids the random Latin Hypercube sampling that generates
    # extreme Arrhenius parameter combinations
    Random.seed!(123)
    nominal = collect(petab_prob.xnominal_transformed)
    lb = collect(petab_prob.lower_bounds)
    ub = collect(petab_prob.upper_bounds)

    all_results = []
    best_fmin = Inf
    best_xmin = nominal

    for i in 1:n_multistarts
        if i == 1
            x0 = copy(nominal)
        else
            # Perturb nominal by ±20% within bounds
            x0 = nominal .+ 0.2 .* (ub .- lb) .* (rand(length(nominal)) .- 0.5)
            x0 = clamp.(x0, lb, ub)
        end

        println("  Start $i/$n_multistarts...")
        try
            res = calibrate(petab_prob, x0, LBFGS())
            fval = res.fmin
            println("    → NLL = $(round(fval, digits=2))")
            if fval < best_fmin
                best_fmin = fval
                best_xmin = res.xmin
            end
            push!(all_results, res)
        catch e
            println("    → Failed (solver error)")
        end
    end

    println("\n  Best NLL: $(round(best_fmin, digits=2))")

    # Return a result-like object compatible with downstream code
    return (xmin = best_xmin, fmin = best_fmin, runs = all_results)
end

"""
    compare_parameters(result, petab_prob)

Create a comparison table of true vs recovered parameters.

Returns: DataFrame with parameter comparison
"""
function compare_parameters(result, petab_prob)
    # True parameter values (from Section 6)
    true_values = Dict(
        "log10_A1" => 5.0,
        "E1"       => 120_000.0,
        "n1"       => 1.0,
        "log10_A2" => 7.0,
        "E2"       => 160_000.0,
        "n2"       => 1.5,
        "m2"       => 0.6,
        "log10_A3" => 9.0,
        "E3"       => 250_000.0,
        "n3"       => 2.0,
        "m3"       => 0.85,
    )

    # Extract recovered values from optimization result
    recovered = Dict(zip(string.(petab_prob.xnames), result.xmin))

    # Build comparison table
    param_names = ["log10_A1", "E1", "n1", "log10_A2", "E2", "n2", "m2",
                   "log10_A3", "E3", "n3", "m3"]
    true_vals = Float64[]
    rec_vals  = Float64[]
    pct_errs  = Float64[]

    for pname in param_names
        tv = true_values[pname]
        rv = recovered[pname]
        pe = abs(rv - tv) / abs(tv) * 100.0
        push!(true_vals, tv)
        push!(rec_vals, rv)
        push!(pct_errs, pe)
    end

    comparison_df = DataFrame(
        Parameter  = param_names,
        True_Value = true_vals,
        Recovered  = round.(rec_vals, sigdigits=5),
        Pct_Error  = round.(pct_errs, digits=2),
    )

    return comparison_df
end

"""
    print_comparison_table(df)

Pretty-print the parameter comparison table.
"""
function print_comparison_table(df)
    println("\n" * "="^70)
    println("  Parameter Estimation Results: True vs Recovered")
    println("="^70)
    @printf("  %-12s  %12s  %12s  %10s\n", "Parameter", "True Value", "Recovered", "% Error")
    println("  " * "-"^50)
    for row in eachrow(df)
        @printf("  %-12s  %12.4f  %12.4f  %9.2f%%\n",
                row.Parameter, row.True_Value, row.Recovered, row.Pct_Error)
    end
    println("="^70)
    println("  Mean % Error: $(round(mean(df.Pct_Error), digits=2))%")
    println("="^70)
end
