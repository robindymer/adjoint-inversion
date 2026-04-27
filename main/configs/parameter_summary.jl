##############################
# Parameter summary function
##############################
function format_parameters(; dt=nothing, hx=nothing, hy=nothing)
    lines = String[]
    push!(lines, "=" ^ 60)
    push!(lines, "SIMULATION PARAMETERS")
    push!(lines, "=" ^ 60)
    
    push!(lines, "\n[Problem Setup]")
    push!(lines, "  Numerical scheme:    $SCHEMA")
    push!(lines, "  SBP order:           $SBP_ORDER")
    push!(lines, "  Precision:           $PRECISION")
    push!(lines, "  Use MMS:             $USE_MMS")
    push!(lines, "  Initial condition:   $INITIAL_CONDITION")
    push!(lines, "  Periodic:            $PERIODIC")
    push!(lines, "  Use CUDA:            $USE_CUDA")
    
    push!(lines, "\n[Grid]")
    push!(lines, "  Grid points:         $mx × $my")
    push!(lines, "  Domain (x):          [$x_l, $x_r]")
    push!(lines, "  Domain (y):          [$y_l, $y_r]")
    if hx !== nothing && hy !== nothing
        push!(lines, "  hx:                  $hx")
        push!(lines, "  hy:                  $hy")
    end
    
    push!(lines, "\n[Boundary Conditions]")
    push!(lines, "  West:                $BC_WEST")
    push!(lines, "  East:                $BC_EAST")
    push!(lines, "  South:               $BC_SOUTH")
    push!(lines, "  North:               $BC_NORTH")
    
    push!(lines, "\n[Material Properties]")
    push!(lines, "  Density (D):         $D_tag")
    push!(lines, "  Lamé λ:              $λ_tag")
    push!(lines, "  Lamé μ:              $μ_tag")
    if USE_MMS
        push!(lines, "  Frequency (ω):       $ω")
    end
    
    push!(lines, "\n[Time Stepping]")
    push!(lines, "  End time (T):        $T")
    if dt !== nothing
        push!(lines, "  Time step (dt):      $dt")
    end
    push!(lines, "  CFL:                 $CFL")
    push!(lines, "  Output interval:     $output_interval")
    
    if PULSE
        push!(lines, "\n[Source/Pulse]")
        push!(lines, "  Enabled:             $PULSE")
        push!(lines, "  Location (X_S):      $X_S")
        push!(lines, "  Center time (T0):    $T0")
        push!(lines, "  Width (SIGMA):       $SIGMA")
        push!(lines, "  Moment order:        $M_ORDER")
        push!(lines, "  Spatial order:       $S_ORDER")
    end
    
    push!(lines, "\n[Output]")
    push!(lines, "  Save VTK:            $SAVE_VTK")
    push!(lines, "  Save GIF:            $SAVE_GIF")
    push!(lines, "  Output path:         $OUTPUT_PATH")
    
    push!(lines, "\n" * "=" ^ 60)
    
    return join(lines, "\n")
end