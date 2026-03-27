# =============================================================================
# CFRP Oxidative Decomposition — Multi-Condition Inverse Modeling
# =============================================================================
# Homework 2: Introduction to Scientific Computing
#
# Complete workflow:
#   Task A — Build Catalyst.jl reaction network
#   Task B — Generate synthetic TGA data (4 experiments)
#   Task C — Set up PEtab inverse problem
#   Task D — Parameter estimation, comparison table, validation plot
#
# Usage:
#   1. Activate the project:  ] activate .
#   2. Install dependencies:  ] instantiate
#   3. Run:                   include("run.jl")
# =============================================================================

using Statistics, CSV

# Include all source modules
include("src/model.jl")
include("src/data_generation.jl")
include("src/inverse_problem.jl")
include("src/visualization.jl")

function main()
    println("="^60)
    println("  CFRP Oxidative Decomposition — Inverse Modeling")
    println("="^60)

    # =========================================================================
    # TASK A: Build the Forward Model
    # =========================================================================
    println("\n📌 Task A: Building Catalyst.jl reaction network...")
    rn = build_cfrp_reaction_system()
    println("  ✓ Reaction system built: $(length(reactions(rn))) reactions, " *
            "$(length(species(rn))) species, $(length(parameters(rn))) parameters")

    # =========================================================================
    # TASK B: Generate Synthetic TGA Data
    # =========================================================================
    println("\n📌 Task B: Generating synthetic TGA data...")
    measurements_df, clean_data = generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)
    println("  ✓ Generated $(nrow(measurements_df)) measurement points across 4 experiments")

    # Save measurement data to CSV
    CSV.write("results/measurements.csv", measurements_df)
    println("  ✓ Measurements saved to results/measurements.csv")

    # Plot all experiments overview
    println("\n  Plotting all experiments...")
    plot_all_experiments(clean_data)

    # =========================================================================
    # TASK C: Set Up the Inverse Problem
    # =========================================================================
    println("\n📌 Task C: Setting up PEtab inverse problem...")
    petab_prob = setup_petab_model(rn, measurements_df)
    println("  ✓ PEtab problem created with $(length(petab_prob.xnames)) parameters to estimate")
    println("  Parameters: ", join(string.(petab_prob.xnames), ", "))

    # =========================================================================
    # TASK D: Parameter Estimation & Validation
    # =========================================================================
    println("\n📌 Task D: Running parameter estimation...")
    result = run_parameter_estimation(petab_prob; n_multistarts=10)

    # --- Comparison Table ---
    comparison_df = compare_parameters(result, petab_prob)
    print_comparison_table(comparison_df)

    # Save comparison table to CSV
    CSV.write("results/parameter_comparison.csv", comparison_df)
    println("  ✓ Comparison table saved to results/parameter_comparison.csv")

    # --- Validation Plot ---
    println("\n  Generating validation plot...")
    plot_validation(rn, result, petab_prob, clean_data)

    # =========================================================================
    # Discussion: Why Experiment 4 is necessary
    # =========================================================================
    println("\n" * "="^70)
    println("  DISCUSSION")
    println("="^70)
    println("""
  Why Experiment 4 (5% O₂) is necessary to determine m₂ and m₃:

  The oxygen reaction orders m₂ and m₃ appear in the rate laws for char
  oxidation and fiber oxidation as r ∝ PO₂^m. If all experiments use the
  same oxygen concentration (21%), then PO₂^m is a constant factor that
  becomes mathematically indistinguishable from the pre-exponential factor
  A — the optimizer can compensate for any m value by adjusting A.
  This creates a parameter identifiability problem where m₂ and m₃ are
  unidentifiable.

  By including Experiment 4 at a different PO₂ (5%), the ratio of reaction
  rates between Exp 2 and Exp 4 at the same temperature depends explicitly
  on m: r(21%)/r(5%) = (0.21/0.05)^m. This additional constraint breaks
  the degeneracy between A and m, making both oxygen reaction orders
  uniquely identifiable from the combined multi-condition dataset.
    """)

    println("="^60)
    println("  ✓ All tasks completed successfully!")
    println("="^60)

    return rn, measurements_df, clean_data, petab_prob, result, comparison_df
end

# Run the full workflow
rn, measurements_df, clean_data, petab_prob, result, comparison_df = main()
