# =============================================================================
# Task A: The Forward Model — CFRP Oxidative Decomposition Reaction Network
# =============================================================================
#
# Three-reaction mechanism for CFRP thermal decomposition:
#   Reaction 1: Matrix Pyrolysis (anaerobic)  — M → ν_char·C + (1-ν_char)·G1
#   Reaction 2: Char Oxidation (aerobic)      — C + O₂ → G2
#   Reaction 3: Fiber Oxidation (aerobic)     — F + O₂ → G2
#
# Temperature modeled as dynamic variable: dT/dt = β (linear heating)
# Rate laws: Arrhenius with power-law dependence on species and O₂ partial pressure
# =============================================================================

using Catalyst, ModelingToolkit

"""
    build_cfrp_reaction_system()

Build the Catalyst.jl ReactionSystem for CFRP oxidative decomposition.

Species: M (Matrix), C (Char), F (Fiber), Temp (Temperature), G1 (Volatiles), G2 (Oxidation gas)
Parameters: log10_A1, E1, n1, nu_char, log10_A2, E2, n2, m2, log10_A3, E3, n3, m3, PO2, beta

Returns a completed ReactionSystem ready for simulation.
"""
function build_cfrp_reaction_system()
    # Independent variable
    t = default_t()

    # Define species
    @species M(t) C(t) F(t) Temp(t) G1(t) G2(t)

    # Define parameters
    # Kinetic parameters (log10 of pre-exponential for numerical stability)
    @parameters log10_A1 E1 n1 nu_char
    @parameters log10_A2 E2 n2 m2
    @parameters log10_A3 E3 n3 m3
    # Experimental condition parameters
    @parameters PO2 beta

    # Gas constant [J/(mol·K)]
    R_gas = 8.314

    # --- Rate expressions ---
    # Use max(species, 0) to prevent NaN from negative^fractional during optimization
    # Reaction 1: Matrix Pyrolysis (anaerobic)
    # r1 = A1 * exp(-E1/(R*T)) * M^n1
    r1 = (10.0^log10_A1) * exp(-E1 / (R_gas * Temp)) * max(M, 0.0)^n1

    # Reaction 2: Char Oxidation (aerobic)
    # r2 = A2 * exp(-E2/(R*T)) * C^n2 * PO2^m2
    r2 = (10.0^log10_A2) * exp(-E2 / (R_gas * Temp)) * max(C, 0.0)^n2 * PO2^m2

    # Reaction 3: Fiber Oxidation (aerobic)
    # r3 = A3 * exp(-E3/(R*T)) * F^n3 * PO2^m3
    r3 = (10.0^log10_A3) * exp(-E3 / (R_gas * Temp)) * max(F, 0.0)^n3 * PO2^m3

    # --- Reactions ---
    rxns = [
        # Temperature ramp: dT/dt = beta  (linear heating program)
        Reaction(beta, nothing, [Temp], nothing, [1]),

        # Reaction 1: M → nu_char*C + (1-nu_char)*G1
        # only_use_rate=true: the rate expression IS the full propensity (not mass-action)
        Reaction(r1, [M], [C, G1], [1], [nu_char, 1 - nu_char]; only_use_rate=true),

        # Reaction 2: C → G2  (with O₂ dependence in rate)
        Reaction(r2, [C], [G2], [1], [1]; only_use_rate=true),

        # Reaction 3: F → G2  (with O₂ dependence in rate)
        Reaction(r3, [F], [G2], [1], [1]; only_use_rate=true),
    ]

    # Build and complete the reaction system
    @named cfrp = ReactionSystem(rxns, t)
    cfrp = complete(cfrp)

    return cfrp
end

# =============================================================================
# True parameters for data generation (Section 6 of assignment)
# =============================================================================
const TRUE_PARAMS = Dict(
    :log10_A1 => 5.0,       # log10(A1) for matrix pyrolysis
    :E1       => 120_000.0,  # Activation energy [J/mol]
    :n1       => 1.0,        # Reaction order for M
    :nu_char  => 0.2,        # Char yield (FIXED, not estimated)
    :log10_A2 => 7.0,        # log10(A2) for char oxidation
    :E2       => 160_000.0,  # Activation energy [J/mol]
    :n2       => 1.5,        # Reaction order for C
    :m2       => 0.6,        # Oxygen reaction order
    :log10_A3 => 9.0,        # log10(A3) for fiber oxidation
    :E3       => 250_000.0,  # Activation energy [J/mol]
    :n3       => 2.0,        # Reaction order for F
    :m3       => 0.85,       # Oxygen reaction order
)

# =============================================================================
# Initial conditions
# =============================================================================
const INITIAL_M    = 0.3    # Matrix mass fraction
const INITIAL_F    = 0.7    # Fiber mass fraction
const INITIAL_C    = 0.0    # Char mass fraction (none initially)
const INITIAL_TEMP = 300.0  # Starting temperature [K]
const INITIAL_G1   = 0.0    # No gas initially
const INITIAL_G2   = 0.0    # No gas initially

# =============================================================================
# Experimental conditions (Section 3)
# Note: β converted from K/min to K/s (kinetics use seconds!)
# =============================================================================
const EXPERIMENTS = Dict(
    "Exp1" => (beta = 2.5 / 60.0,  PO2 = 0.21),  # 2.5 K/min, 21% O₂
    "Exp2" => (beta = 5.0 / 60.0,  PO2 = 0.21),  # 5.0 K/min, 21% O₂
    "Exp3" => (beta = 10.0 / 60.0, PO2 = 0.21),  # 10.0 K/min, 21% O₂
    "Exp4" => (beta = 5.0 / 60.0,  PO2 = 0.05),  # 5.0 K/min, 5% O₂
)

# Final temperature for TGA simulation [K]
const T_FINAL = 1200.0
