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

const USE_CUDA = false
const PRECISION = Float64  # Float32 or Float64

const TEST_SYMMETRY = false
const DEBUG_CUDA = false

const USE_MARMOUSI = true
const DATA_DIR = joinpath(@__DIR__, "../../../data/marmousi2")
const DOWNSAMPLE_FACTOR = 5

# --- Boundary conditions ---
const BC_WEST = :non_reflecting      # :traction_free or :non_reflecting
const BC_EAST = :non_reflecting
const BC_SOUTH = :non_reflecting
const BC_NORTH = :traction_free

# --- Pulse parameters ---
const PULSE = true
const X_S = [1.25*1000, 1.25*15]
const T0 = 1.0
const SIGMA = 0.1149
const M_ORDER = 4
const S_ORDER = 4
const SCALE_GAUSSIAN = true

# --- Grid parameters ---
const mx = 13601 ÷ DOWNSAMPLE_FACTOR                     # grid points in x
const my = 2801 ÷ DOWNSAMPLE_FACTOR                     # grid points in y
const x_l = 0.0
const x_r = 13600*1.25
const y_l = 0.0
const y_r = 2800*1.25

# --- Time stepping ---
const CFL = 0.5
const T = 15.0
const output_interval = 50

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
const SAVE_VTK = true
const SAVE_GIF = false
const ASSETS_PATH = joinpath(@__DIR__, "../../../assets")
const DESC_STRING = "CUDA_$(USE_CUDA)_M-S=$(M_ORDER)-$(S_ORDER)_mumarmousi_t0_$(replace(string(T0), "." => "_"))_SIGMA_$(replace(string(SIGMA), "." => "_"))_SCHEMA_$(SCHEMA)_$(PERIODIC ? "periodic" : "non_periodic")_IC_$(INITIAL_CONDITION)_order_$(SBP_ORDER)_mxmy_$(mx)x$(my)"
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
