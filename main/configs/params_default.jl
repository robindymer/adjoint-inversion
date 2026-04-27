##############################
# Simulation Configuration
##############################

module Params
include("parameter_summary.jl")

# --- Problem setup ---
const USE_MMS = false              # manufactured solution toggle
const INITIAL_CONDITION = :gaussian # :gaussian or :curl_free, NOT relevant if USE_MMS=true
const PERIODIC = false          # periodicity toggle
const SCHEMA = :mixed_upwind_λ            # :narrow, :wide, :mixed, :upwind, :mixed_upwind_μ, :mixed_upwind_λ and :staggered

const USE_CUDA = false
const PRECISION = Float64  # Float32 or Float64

const TEST_SYMMETRY = false
const DEBUG_CUDA = false

const USE_MARMOUSI = false
const DATA_DIR = joinpath(@__DIR__, "../../../data/marmousi2")
const DOWNSAMPLE_FACTOR = 4

# --- Boundary conditions ---
const BC_WEST = :traction_free      # :traction_free or :non_reflecting
const BC_EAST = :traction_free
const BC_SOUTH = :traction_free
const BC_NORTH = :traction_free
# const BC_WEST = :non_reflecting
# const BC_EAST = :non_reflecting
# const BC_SOUTH = :non_reflecting
# const BC_NORTH = :non_reflecting

# --- Pulse parameters ---
const PULSE = false
const X_S = [0.0, 0.0]
const T0 = 1.4
const SIGMA = 0.1149
const M_ORDER = 4
const S_ORDER = 4

# --- Grid parameters ---
const x_l = -1.0
const x_r = 1.0
const y_l = -1.0
const y_r = 1.0
const mx = 100                     # grid points in x
const my = 100                     # grid points in y

# --- Time stepping ---
const CFL = 0.8
const T = 1.0
const output_interval = 10

# --- Material parameters ---
const ω = 1.0                     # frequency for MMS
const λ_tag = "trigonometric"
const μ_tag = "trigonometric"
const D_tag = "trigonometric"

# --- Boundary conditions --- NOT YET SUPPORTED
# const BC_type_x = :traction_free           # :traction_free, :fixed, :periodic
# const BC_type_y = :traction_free

# --- Stencil options ---
const SBP_ORDER = 4                   # SBP order

# --- Output options ---
const SAVE_VTK = false
const SAVE_GIF = true
const ASSETS_PATH = joinpath(@__DIR__, "../assets")
const OUTPUT_PATH = joinpath(@__DIR__, "../../../output/results_SCHEMA_$(SCHEMA)_IC_$(INITIAL_CONDITION)_order_$(SBP_ORDER)")
const OUTPUT_NAME = "$OUTPUT_PATH/simulation_$(SCHEMA)_$(INITIAL_CONDITION)_order_$(SBP_ORDER)"
const VERBOSE = true

# --- Convergence study options ---
const CONVERGENCE_STUDY = false
const m_values = [20, 40, 80]

##############################
# Print summary if verbose
##############################
if VERBOSE
    println(format_parameters())
end

end # module
