using JLD2
using Printf

println("="^60)
println("GPU vs CPU Derivative Computation Comparison")
println("="^60)

# Load data
println("\nLoading data files...")
# gpu_inputs = load("cuda_kernel_inputs_GPU.jld2")
gpu_outputs = load("cuda_kernel_outputs_GPU.jld2")
# cpu_inputs = load("cuda_kernel_inputs_CPU.jld2")
cpu_outputs = load("cuda_kernel_outputs_CPU.jld2")

println("✓ Files loaded successfully")

# # Compare inputs first (sanity check)
# println("\n" * "="^60)
# println("INPUTS COMPARISON (Sanity Check)")
# println("="^60)

# input_vars = ["u1", "u2", "mx", "my", "Ω_radius", "λ", "μ", "D1x", "D1y", "BW_D1"]

# all_inputs_match = true
# for var in input_vars
#     if haskey(gpu_inputs, var) && haskey(cpu_inputs, var)
#         gpu_val = gpu_inputs[var]
#         cpu_val = cpu_inputs[var]
        
#         if isa(gpu_val, Array) && isa(cpu_val, Array)
#             match = gpu_val ≈ cpu_val
#             if match
#                 println("✓ $var: MATCH")
#             else
#                 max_diff = maximum(abs.(gpu_val .- cpu_val))
#                 max_rel_diff = maximum(abs.((gpu_val .- cpu_val) ./ (abs.(cpu_val) .+ 1e-10)))
#                 println("✗ $var: MISMATCH")
#                 println("  Max absolute diff: $max_diff")
#                 println("  Max relative diff: $max_rel_diff")
#                 global all_inputs_match = false
                
#                 # Show some sample differences
#                 diff = abs.(gpu_val .- cpu_val)
#                 worst_idx = argmax(diff)
#                 println("  Worst mismatch at index $worst_idx:")
#                 println("    GPU: $(gpu_val[worst_idx])")
#                 println("    CPU: $(cpu_val[worst_idx])")
#             end
#         else
#             match = gpu_val == cpu_val
#             println(match ? "✓ $var: $gpu_val" : "✗ $var: GPU=$gpu_val, CPU=$cpu_val")
#             global all_inputs_match = all_inputs_match && match
#         end
#     end
# end

# if all_inputs_match
#     println("\n✓ All inputs match - good starting point!")
# else
#     println("\n⚠ WARNING: Inputs don't match! Results may not be comparable.")
# end

# Compare outputs
println("\n" * "="^60)
println("OUTPUTS COMPARISON")
println("="^60)

# output_vars = ["u1_x", "u1_y", "u2_x", "u2_y", "u2_yλx", "u2_xμy", "u1_xλy", "u1_yμx"]
output_vars = ["u1_x", "u1_y", "u2_x", "u2_y", "sigma_11_x", "sigma_11_y", "sigma_12_x", "sigma_12_y", "sigma_22_x", "sigma_22_y"]

all_outputs_match = true
for var in output_vars
    if haskey(gpu_outputs, var) && haskey(cpu_outputs, var)
        gpu_val = gpu_outputs[var]
        cpu_val = cpu_outputs[var]
        
        println("\n" * "-"^60)
        println("Variable: $var")
        println("-"^60)
        
        # Check dimensions
        if size(gpu_val) != size(cpu_val)
            println("✗ SIZE MISMATCH: GPU=$(size(gpu_val)), CPU=$(size(cpu_val))")
            global all_outputs_match = false
            continue
        end
        
        # Calculate statistics
        diff = abs.(gpu_val .- cpu_val)
        max_abs_diff = maximum(diff)
        mean_abs_diff = sum(diff) / length(diff)
        
        # Relative error (avoid division by zero)
        rel_diff = abs.((gpu_val .- cpu_val) ./ (abs.(cpu_val) .+ 1e-10))
        max_rel_diff = maximum(rel_diff)
        mean_rel_diff = sum(rel_diff) / length(rel_diff)
        
        # Count zeros
        gpu_zeros = count(x -> abs(x) < 1e-10, gpu_val)
        cpu_zeros = count(x -> abs(x) < 1e-10, cpu_val)
        
        # Count non-zeros that differ significantly
        significant_diffs = count(diff .> 1e-6)
        
        match = max_abs_diff < 1e-6
        
        println("  Size: $(size(gpu_val))")
        println("  GPU zeros: $gpu_zeros / $(length(gpu_val))")
        println("  CPU zeros: $cpu_zeros / $(length(cpu_val))")
        println("  Max absolute diff: $(@sprintf("%.6e", max_abs_diff))")
        println("  Mean absolute diff: $(@sprintf("%.6e", mean_abs_diff))")
        println("  Max relative diff: $(@sprintf("%.6e", max_rel_diff))")
        println("  Mean relative diff: $(@sprintf("%.6e", mean_rel_diff))")
        println("  Significant diffs (>1e-6): $significant_diffs")
        
        if match
            println("  ✓ MATCH (within tolerance)")
        else
            println("  ✗ MISMATCH")
            global all_outputs_match = false
            
            # Find and show worst mismatches
            worst_indices = sortperm(vec(diff), rev=true)[1:min(5, length(diff))]
            println("\n  Top 5 worst mismatches:")
            for (i, idx) in enumerate(worst_indices)
                println("    $i. Index $idx:")
                println("       GPU: $(@sprintf("%.10e", gpu_val[idx]))")
                println("       CPU: $(@sprintf("%.10e", cpu_val[idx]))")
                println("       Diff: $(@sprintf("%.10e", diff[idx]))")
            end
            
            # Show distribution of differences
            println("\n  Difference distribution:")
            bins = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0, Inf]
            for i in 1:length(bins)-1
                count_in_bin = count(x -> bins[i] <= x < bins[i+1], diff)
                pct = 100 * count_in_bin / length(diff)
                println("    [$(@sprintf("%.0e", bins[i])), $(@sprintf("%.0e", bins[i+1]))): $count_in_bin ($(@sprintf("%.1f", pct))%)")
            end
        end
    else
        println("\n✗ $var: Missing in one of the output files")
        global all_outputs_match = false
    end
end

# Summary
println("\n" * "="^60)
println("SUMMARY")
println("="^60)

if all_outputs_match
    println("✓ All outputs match within tolerance!")
else
    println("✗ Outputs differ between GPU and CPU implementations")
    println("\nPossible issues to check:")
    println("  1. Indexing bugs in CUDA kernels")
    println("  2. Race conditions in atomic operations")
    println("  3. Boundary treatment differences")
    println("  4. Floating point precision issues (Float32 vs Float64)")
end

println("\n" * "="^60)
