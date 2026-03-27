# CFRP Oxidative Decomposition — Multi-Condition Inverse Modeling

Homework 2 for Introduction to Scientific Computing (Dr. Hayri Sezer) class.

## Problem

Multi-stage kinetic modeling of Carbon Fiber Reinforced Polymer (CFRP) thermal decomposition:

| Reaction | Description | Rate Law |
|----------|-------------|----------|
| 1. Matrix Pyrolysis | M → ν_char·C + (1-ν_char)·G1 | r₁ = A₁·exp(-E₁/RT)·[M]^n₁ |
| 2. Char Oxidation | C + O₂ → G₂ | r₂ = A₂·exp(-E₂/RT)·[C]^n₂·PO₂^m₂ |
| 3. Fiber Oxidation | F + O₂ → G₂ | r₃ = A₃·exp(-E₃/RT)·[F]^n₃·PO₂^m₃ |

## Project Structure

```
cfrp_kinetics/
├── Project.toml            # Julia dependencies
├── README.md
├── run.jl                  # Main entry point (runs full workflow)
├── src/
│   ├── model.jl            # Task A: Catalyst.jl reaction network
│   ├── data_generation.jl  # Task B: Synthetic TGA data generation
│   ├── inverse_problem.jl  # Tasks C & D: PEtab setup + optimization
│   └── visualization.jl    # Validation plots
└── results/                # Output directory (plots, CSVs)
```

## How to Run

```julia
# 1. Open Julia in the project directory
cd("path/to/cfrp_kinetics")

# 2. Activate and install dependencies
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# 3. Run the full workflow
include("run.jl")
```

## Experiments

| Exp | β (K/min) | PO₂ | Purpose |
|-----|-----------|------|---------|
| 1 | 2.5 | 0.21 | Kinetic Baseline |
| 2 | 5.0 | 0.21 | Kinetic Baseline |
| 3 | 10.0 | 0.21 | Determine Activation Energy |
| 4 | 5.0 | 0.05 | Determine Oxygen Order (m) |

## Dependencies

- Catalyst.jl — Reaction network definition
- PEtab.jl — Parameter estimation (inverse problem)
- OrdinaryDiffEq.jl — ODE solver
- Plots.jl — Visualization
- DataFrames.jl, CSV.jl — Data handling
