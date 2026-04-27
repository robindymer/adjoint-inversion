##############################
# Simulation Configuration
##############################
module Params
include("parameter_summary.jl")

# --- Problem setup ---
const USE_MMS = false              # manufactured solution toggle
const INITIAL_CONDITION = :zero # :gaussian, :curl_free or :zero, NOT relevant if USE_MMS=true
const PERIODIC = false          # periodicity toggle
const SCHEMA = length(ARGS) > 1 ? Symbol(ARGS[2]) : :narrow           # :narrow, :wide, :mixed, :upwind, :mixed_upwind_μ, :mixed_upwind_λ and :staggered

const USE_CUDA = true
const PRECISION = Float32  # Float32 or Float64

const TEST_SYMMETRY = false
const DEBUG_CUDA = false

const USE_MARMOUSI = false
const DATA_DIR = joinpath(@__DIR__, "../../../data/marmousi2")
const DOWNSAMPLE_FACTOR = 4

# --- Boundary conditions ---
const BC_WEST = :non_reflecting      # :traction_free or :non_reflecting
const BC_EAST = :non_reflecting
const BC_SOUTH = :non_reflecting
const BC_NORTH = :non_reflecting

# --- Pulse parameters ---
const PULSE = true
const X_S = [0.0, 0.0]
const T0 = 1.5
const SIGMA = 0.2
const M_ORDER = 4
const S_ORDER = 4

# --- Grid parameters ---
const x_l = -5.0
const x_r = 5.0
const y_l = -5.0
const y_r = 5.0
const mx = 100                     # grid points in x
const my = 100                     # grid points in y

# --- Time stepping ---
const CFL = 0.4
const T = 12.0
const output_interval = 100

# --- Material parameters ---
const ω = 1.0                     # frequency for MMS
const λ_tag = "constant"
const μ_tag = "constant_005"
const D_tag = "constant"

# --- Boundary conditions --- NOT YET SUPPORTED
# const BC_type_x = :traction_free           # :traction_free, :fixed, :periodic
# const BC_type_y = :traction_free

# --- Stencil options ---
const SBP_ORDER = 4

# --- Output options ---
const SAVE_VTK = false
const SAVE_GIF = true
const ASSETS_PATH = joinpath(@__DIR__, "../../../assets")
const DESC_STRING = "CUDA_M-S=$(M_ORDER)-$(S_ORDER)_mu$(μ_tag)_t0_$(replace(string(T0), "." => "_"))_SIGMA_$(replace(string(SIGMA), "." => "_"))_SCHEMA_$(SCHEMA)_$(PERIODIC ? "periodic" : "non_periodic")_IC_$(INITIAL_CONDITION)_order_$(SBP_ORDER)_mxmy_$(mx)x$(my)"
const OUTPUT_PATH = joinpath(@__DIR__, "../../../output/results_$(DESC_STRING)")
const OUTPUT_NAME = "$OUTPUT_PATH/simulation_$(DESC_STRING)"
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
