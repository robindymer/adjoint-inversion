# Gaussian initial condition or Curl free initial condition
function initialcondition(x, y; a=0.2)
    if Params.INITIAL_CONDITION == :curl_free
        φ = exp(-((x^2 + y^2) / a^2))
        ux = -2x / a^2 * φ
        uy = -2y / a^2 * φ
        return (ux, uy)
    elseif Params.INITIAL_CONDITION == :gaussian
        ux = exp(-((x).^2 + (y).^2) / (a^2))
        uy = exp(-((x).^2 + (y).^2) / (a^2))
        return (ux, uy)
    elseif Params.INITIAL_CONDITION == :zero
        return (0.0, 0.0)
    else
        error("Unknown initial condition type")
    end
end

function sourceTime(t::Float64, t0::Float64, sigma::Float64)
    # Gaussian pulse
    if Params.SCALE_GAUSSIAN
    	return 1/(sigma*sqrt(2*pi))*exp(-(t-t0)^2/(2*sigma^2))
    else
	    return exp(-(t-t0)^2/(2*sigma^2))
    end
    # Ricker wavelet
   # return 2/(sqrt(3*sigma)*pi^(1/4))*(1-((t - t0)^2)/(sigma^2))*exp(-((t - t0)^2)/(2*sigma^2))
end
