### CUDA kernels in this file

# TODO: Fix code duplication, try e.g. if Params.SCHEMA and hope that warp divergence is not too bad
function apply_narrow_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, D1x, D1y, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    # -- Good when the dataset size exceeds the total number of threads in the grid --
    # stridei = blockDim().x * gridDim().x
    # stridej = blockDim().y * gridDim().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1
    # D2x_boundary only sees x direction, D2y_boundary only y direction
    hx2inv = 1.0f0/(hx*hx)
        D2x_lambda = -@SVector [
        hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
        hx2inv*(-lambda[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
        hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + lambda[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
        hx2inv*(-lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
        hx2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    ]
    D2x_mu = -@SVector [
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
        hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + mu[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
        hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
        hx2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + mu[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    ]
    hy2inv = 1.0f0/(hy*hy)
    D2y_mu = -@SVector [
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
        hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + mu[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
        hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
        hy2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + mu[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    ]
    D2y_lambda = -@SVector [
        hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
        hy2inv*(-lambda[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
        hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + lambda[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
        hy2inv*(-lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
        hy2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    ]

    interior_offset = Omega_radius ÷ 2

    # for l = i:stridei:mx
    #     for m = j:stridej:my
    #     end
    # end

    # -- NOTE: While loops faster, but more bloated code --
    # k1 = 1
    # while k1 <= Omega_radius
    for k1=1:Omega_radius
        # k2 = 1
        wrapped_x = wrap_x(i - interior_offset - 1 + k1)
        wrapped_y = wrap_y(j - interior_offset - 1 + k1)
        # -- D1 --
        u1_x[idx(i, j, mx, my)] += D1x[k1] * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k1] * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k1] * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k1] * u2[idx(i, wrapped_y, mx, my)]

        # -- D2, narrow --
        u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda[k1] * u1[idx(wrapped_x, j, mx, my)]
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k1] * u1[idx(wrapped_x, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k1] * u2[idx(wrapped_x, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k1] * u1[idx(i, wrapped_y, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k1] * u2[idx(i, wrapped_y, mx, my)]
        u2_ylambday[idx(i, j, mx, my)] += D2y_lambda[k1] * u2[idx(i, wrapped_y, mx, my)]

        # -- D2, mixed --
        # while k2 <= Omega_radius
        for k2=1:Omega_radius
            # Periodic wrapping
            wrapped_x = wrap_x(i - interior_offset - 1 + k1)
            wrapped_y = wrap_y(j - interior_offset - 1 + k2)
            u2_ylambdax[idx(i, j, mx, my)] += D1x[k1] * lambda[idx(wrapped_x, j, mx, my)] * D1y[k2] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y[k2] * mu[idx(i, wrapped_y, mx, my)] * D1x[k1] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y[k2] * lambda[idx(i, wrapped_y, mx, my)] * D1x[k1] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x[k1] * mu[idx(wrapped_x, j, mx, my)] * D1y[k2] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            # k2 += 1
        end

        # k1 += 1
    end

    return nothing
end

function apply_wide_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, D1x, D1y, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    # -- Good when the dataset size exceeds the total number of threads in the grid --
    # stridei = blockDim().x * gridDim().x
    # stridej = blockDim().y * gridDim().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1

    interior_offset = Omega_radius ÷ 2

    # -- D1 --
    for k=1:Omega_radius
        wrapped_x = wrap_x(i - interior_offset - 1 + k)
        wrapped_y = wrap_y(j - interior_offset - 1 + k)
        # -- D1 --
        u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, wrapped_y, mx, my)]
    end
    # -- D2, wide (apply D1 twice with material coefficients in between) --
    # x-direction: u1_xlambdax = ∂x(lambda ∂x u1), u1_xmux = ∂x(mu ∂x u1), u2_xmux = ∂x(mu ∂x u2)
    for kx=1:Omega_radius
        wrapped_kx = wrap_x(i - interior_offset - 1 + kx)

        # Compute inner derivative: ∂x u at position wrapped_kx
        inner_u1_x = 0.0
        inner_u2_x = 0.0
        @inbounds for m=1:Omega_radius
            wrapped_m = wrap_x(wrapped_kx - interior_offset - 1 + m)
            inner_u1_x += D1x[m] * u1[idx(wrapped_m, j, mx, my)]
            inner_u2_x += D1x[m] * u2[idx(wrapped_m, j, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = D1x[kx]
        lambdakx = lambda[idx(wrapped_kx, j, mx, my)]
        mukx = mu[idx(wrapped_kx, j, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += w * lambdakx * inner_u1_x
        u1_xmux[idx(i, j, mx, my)] += w * mukx * inner_u1_x
        u2_xmux[idx(i, j, mx, my)] += w * mukx * inner_u2_x
    end

    # y-direction: u1_ymuy = ∂y(mu ∂y u1), u2_ymuy = ∂y(mu ∂y u2), u2_ylambday = ∂y(lambda ∂y u2)
    for ky=1:Omega_radius
        wrapped_ky = wrap_y(j - interior_offset - 1 + ky)

        # Compute inner derivative: ∂y u at position wrapped_ky
        inner_u1_y = 0.0
        inner_u2_y = 0.0
        @inbounds for n=1:Omega_radius
            wrapped_n = wrap_y(wrapped_ky - interior_offset - 1 + n)
            inner_u1_y += D1y[n] * u1[idx(i, wrapped_n, mx, my)]
            inner_u2_y += D1y[n] * u2[idx(i, wrapped_n, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = D1y[ky]
        lambdaky = lambda[idx(i, wrapped_ky, mx, my)]
        muky = mu[idx(i, wrapped_ky, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += w * muky * inner_u1_y
        u2_ymuy[idx(i, j, mx, my)] += w * muky * inner_u2_y
        u2_ylambday[idx(i, j, mx, my)] += w * lambdaky * inner_u2_y
    end

    for kx=1:Omega_radius
        for ky=1:Omega_radius
            # Periodic wrapping
            wrapped_x = wrap_x(i - interior_offset - 1 + kx)
            wrapped_y = wrap_y(j - interior_offset - 1 + ky)
            u2_ylambdax[idx(i, j, mx, my)] += D1x[kx] * lambda[idx(wrapped_x, j, mx, my)] * D1y[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, wrapped_y, mx, my)] * D1x[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, wrapped_y, mx, my)] * D1x[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x[kx] * mu[idx(wrapped_x, j, mx, my)] * D1y[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
        end
    end

    return nothing
end

function apply_mixed_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, D1x, D1y, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    # -- Good when the dataset size exceeds the total number of threads in the grid --
    # stridei = blockDim().x * gridDim().x
    # stridej = blockDim().y * gridDim().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1

    # D2x_boundary only sees x direction, D2y_boundary only y direction
    hx2inv = 1.0f0/(hx*hx)
    D2x_mu = -@SVector [
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/0.8e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/0.6e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/2.4e1 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2e0 + mu[idx(wrap_x(i+2), j, mx, my)]/2.4e1),
        hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+2), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1 + mu[idx(wrap_x(i+2), j, mx, my)]/0.8e1)
    ]
    hy2inv = 1.0f0/(hy*hy)
    D2y_mu = -@SVector [
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/0.8e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/0.6e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/2.4e1 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2e0 + mu[idx(i, wrap_y(j+2), mx, my)]/2.4e1),
        hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+2), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1 + mu[idx(i, wrap_y(j+2), mx, my)]/0.8e1)
    ]

    interior_offset = Omega_radius ÷ 2

    # -- D2, mixed i.e. wide for lambda and narrow for mu --
    for k=1:Omega_radius
        wrapped_x = wrap_x(i - interior_offset - 1 + k)
        wrapped_y = wrap_y(j - interior_offset - 1 + k)
        # -- D1 --
        u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, wrapped_y, mx, my)]

        # -- D2, narrow --
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(wrapped_x, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(wrapped_x, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, wrapped_y, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, wrapped_y, mx, my)]
    end
    for kx=1:Omega_radius
        wrapped_kx = wrap_x(i - interior_offset - 1 + kx)

        # Compute inner derivative: ∂x u at position wrapped_kx
        inner_u1_x = 0.0
        @inbounds for m=1:Omega_radius
            wrapped_m = wrap_x(wrapped_kx - interior_offset - 1 + m)
            inner_u1_x += D1x[m] * u1[idx(wrapped_m, j, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = D1x[kx]
        lambdakx = lambda[idx(wrapped_kx, j, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += w * lambdakx * inner_u1_x
    end

    for ky=1:Omega_radius
        wrapped_ky = wrap_y(j - interior_offset - 1 + ky)

        # Compute inner derivative: ∂y u at position wrapped_ky
        inner_u2_y = 0.0
        @inbounds for n=1:Omega_radius
            wrapped_n = wrap_y(wrapped_ky - interior_offset - 1 + n)
            inner_u2_y += D1y[n] * u2[idx(i, wrapped_n, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = D1y[ky]
        lambdaky = lambda[idx(i, wrapped_ky, mx, my)]

        u2_ylambday[idx(i, j, mx, my)] += w * lambdaky * inner_u2_y
    end

    for kx=1:Omega_radius
        for ky=1:Omega_radius
            # Periodic wrapping
            wrapped_x = wrap_x(i - interior_offset - 1 + kx)
            wrapped_y = wrap_y(j - interior_offset - 1 + ky)
            u2_ylambdax[idx(i, j, mx, my)] += D1x[kx] * lambda[idx(wrapped_x, j, mx, my)] * D1y[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, wrapped_y, mx, my)] * D1x[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, wrapped_y, mx, my)] * D1x[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x[kx] * mu[idx(wrapped_x, j, mx, my)] * D1y[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
        end
    end

    return nothing
end

function apply_upwind_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, Dpx, Dpy, Dmx, Dmy, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1

    interior_offset = Omega_radius ÷ 2

    # elseif SCHEMA == :upwind
    # -- D1 with centered stencil (average of Dp and Dm) --
    for k=1:Omega_radius
        wrapped_x = wrap_x(i - interior_offset - 1 + k)
        wrapped_y = wrap_y(j - interior_offset - 1 + k)
        # Use average of upwind operators to get centered derivative
        D1x_centered = 0.5 * (Dpx[k] + Dmx[k])
        D1y_centered = 0.5 * (Dpy[k] + Dmy[k])

        u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
    end

    # -- D2, wide: Use Dm(c * Dp(u)) for SBP property --
    # x-direction: u1_xlambdax = ∂x(lambda ∂x u1), u1_xmux = ∂x(mu ∂x u1), u2_xmux = ∂x(mu ∂x u2)
    for kx=1:Omega_radius
        wrapped_kx = wrap_x(i - interior_offset - 1 + kx)

        # Compute inner derivative: ∂x u at position wrapped_kx
        inner_u1_x = 0.0
        inner_u2_x = 0.0
        @inbounds for m=1:Omega_radius
            wrapped_m = wrap_x(wrapped_kx - interior_offset - 1 + m)
            inner_u1_x += Dmx[m] * u1[idx(wrapped_m, j, mx, my)]
            inner_u2_x += Dmx[m] * u2[idx(wrapped_m, j, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = Dpx[kx]
        lambdakx = lambda[idx(wrapped_kx, j, mx, my)]
        mukx = mu[idx(wrapped_kx, j, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += w * lambdakx * inner_u1_x
        u1_xmux[idx(i, j, mx, my)] += w * mukx * inner_u1_x
        u2_xmux[idx(i, j, mx, my)] += w * mukx * inner_u2_x
    end

    # y-direction: u1_ymuy = ∂y(mu ∂y u1), u2_ymuy = ∂y(mu ∂y u2), u2_ylambday = ∂y(lambda ∂y u2)
    for ky=1:Omega_radius
        wrapped_ky = wrap_y(j - interior_offset - 1 + ky)

        # Compute inner derivative: ∂y u at position wrapped_ky
        inner_u1_y = 0.0
        inner_u2_y = 0.0
        @inbounds for n=1:Omega_radius
            wrapped_n = wrap_y(wrapped_ky - interior_offset - 1 + n)
            inner_u1_y += Dmy[n] * u1[idx(i, wrapped_n, mx, my)]
            inner_u2_y += Dmy[n] * u2[idx(i, wrapped_n, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = Dpy[ky]
        lambdaky = lambda[idx(i, wrapped_ky, mx, my)]
        muky = mu[idx(i, wrapped_ky, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += w * muky * inner_u1_y
        u2_ymuy[idx(i, j, mx, my)] += w * muky * inner_u2_y
        u2_ylambday[idx(i, j, mx, my)] += w * lambdaky * inner_u2_y
    end

    for kx=1:Omega_radius
        for ky=1:Omega_radius
            # Periodic wrapping
            wrapped_x = wrap_x(i - interior_offset - 1 + kx)
            wrapped_y = wrap_y(j - interior_offset - 1 + ky)

            u2_ylambdax[idx(i, j, mx, my)] += Dpx[kx] * lambda[idx(wrapped_x, j, mx, my)] * Dmy[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += Dpy[ky] * mu[idx(i, wrapped_y, mx, my)] * Dmx[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += Dpy[ky] * lambda[idx(i, wrapped_y, mx, my)] * Dmx[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += Dpx[kx] * mu[idx(wrapped_x, j, mx, my)] * Dmy[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
        end
    end

    return nothing
end

function apply_mixed_mu_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, narrow_Omega_radius, Dpx, Dpy, Dmx, Dmy, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1

    hx2inv = 1.0f0/(hx*hx)
    D2x_mu = -@SVector [
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/0.8e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/0.6e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/2.4e1 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2e0 + mu[idx(wrap_x(i+2), j, mx, my)]/2.4e1),
        hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+2), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1 + mu[idx(wrap_x(i+2), j, mx, my)]/0.8e1)
    ]
    hy2inv = 1.0f0/(hy*hy)
    D2y_mu = -@SVector [
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/0.8e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/0.6e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/2.4e1 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2e0 + mu[idx(i, wrap_y(j+2), mx, my)]/2.4e1),
        hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+2), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1 + mu[idx(i, wrap_y(j+2), mx, my)]/0.8e1)
    ]

    interior_offset = Omega_radius ÷ 2

    # elseif SCHEMA == :mixed_upwind_mu
    narrow_Omega_radius = narrow_Omega_radius
    narrow_interior_offset = narrow_Omega_radius ÷ 2
    # -- D2, mixed i.e. wide for lambda and narrow for mu --
    for k=1:Omega_radius
        wrapped_x = wrap_x(i - interior_offset - 1 + k)
        wrapped_y = wrap_y(j - interior_offset - 1 + k)
        # -- D1 --
        D1x_centered = 0.5 * (Dpx[k] + Dmx[k])
        D1y_centered = 0.5 * (Dpy[k] + Dmy[k])
        # NOTE: Can do it like this instead
        # D1x = [1/12, -2/3, 0, 2/3, -1/12]*1/grid.hx
        # D1y = [1/12, -2/3, 0, 2/3, -1/12]*1/grid.hy
        # D1x_centered = D1x[k]
        # D1y_centered = D1y[k]

        u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
    end
    for k=1:narrow_Omega_radius
        wrapped_x = wrap_x(i - narrow_interior_offset - 1 + k)
        wrapped_y = wrap_y(j - narrow_interior_offset - 1 + k)

        # -- D2, narrow --
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(wrapped_x, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(wrapped_x, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, wrapped_y, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, wrapped_y, mx, my)]
    end

    # Wide derivatives with centered stencils
    # Compute centered stencils on-the-fly to avoid array allocation in GPU kernel
    for kx=1:Omega_radius
        wrapped_kx = wrap_x(i - interior_offset - 1 + kx)

        # Compute inner derivative: ∂x u at position wrapped_kx
        inner_u1_x = 0.0
        @inbounds for m=1:Omega_radius
            wrapped_m = wrap_x(wrapped_kx - interior_offset - 1 + m)
            D1x_centered_m = 0.5 * (Dpx[m] + Dmx[m])
            inner_u1_x += D1x_centered_m * u1[idx(wrapped_m, j, mx, my)]
        end

        # Apply outer derivative with material coefficients
        D1x_centered_kx = 0.5 * (Dpx[kx] + Dmx[kx])
        lambdakx = lambda[idx(wrapped_kx, j, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += D1x_centered_kx * lambdakx * inner_u1_x
    end

    for ky=1:Omega_radius
        wrapped_ky = wrap_y(j - interior_offset - 1 + ky)

        # Compute inner derivative: ∂y u at position wrapped_ky
        inner_u2_y = 0.0
        @inbounds for n=1:Omega_radius
            wrapped_n = wrap_y(wrapped_ky - interior_offset - 1 + n)
            D1y_centered_n = 0.5 * (Dpy[n] + Dmy[n])
            inner_u2_y += D1y_centered_n * u2[idx(i, wrapped_n, mx, my)]
        end

        # Apply outer derivative with material coefficients
        D1y_centered_ky = 0.5 * (Dpy[ky] + Dmy[ky])
        lambdaky = lambda[idx(i, wrapped_ky, mx, my)]

        u2_ylambday[idx(i, j, mx, my)] += D1y_centered_ky * lambdaky * inner_u2_y
    end

    for kx=1:Omega_radius
        for ky=1:Omega_radius
            wrapped_x = wrap_x(i - interior_offset - 1 + kx)
            wrapped_y = wrap_y(j - interior_offset - 1 + ky)

            # Centered stencils from averaging upwind operators
            D1x_centered = 0.5 * (Dpx[kx] + Dmx[kx])
            D1y_centered = 0.5 * (Dpy[ky] + Dmy[ky])

            u2_ylambdax[idx(i, j, mx, my)] += D1x_centered * lambda[idx(wrapped_x, j, mx, my)] * D1y_centered * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += Dpy[ky] * mu[idx(i, wrapped_y, mx, my)] * Dmx[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y_centered * lambda[idx(i, wrapped_y, mx, my)] * D1x_centered * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += Dpx[kx] * mu[idx(wrapped_x, j, mx, my)] * Dmy[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
        end
    end

    return nothing
end

function apply_mixed_lambda_kernel!(u1_x, u1_y, u2_x, u2_y, u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                            u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                            u1, u2, mx, my, hx, hy, Omega_radius, narrow_Omega_radius, Dpx, Dpy, Dmx, Dmy, lambda, mu)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    @inline wrap_y(idx) = (idx+my-1) % my + 1

    hx2inv = 1.0f0/(hx*hx)
    D2x_mu = -@SVector [
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/0.8e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/0.6e1 - mu[idx(wrap_x(i-1), j, mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/2.4e1 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2e0 + mu[idx(wrap_x(i+2), j, mx, my)]/2.4e1),
        hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.2e1 - mu[idx(wrap_x(i+2), j, mx, my)]/0.6e1),
        hx2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(wrap_x(i+1), j, mx, my)]/0.6e1 + mu[idx(wrap_x(i+2), j, mx, my)]/0.8e1)
    ]
    hy2inv = 1.0f0/(hy*hy)
    D2y_mu = -@SVector [
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/0.8e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 + mu[idx(i, j, mx, my)]/0.8e1),
        hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/0.6e1 - mu[idx(i, wrap_y(j-1), mx, my)]/0.2e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/2.4e1 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2e0 + mu[idx(i, j, mx, my)]*0.3/0.4 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2e0 + mu[idx(i, wrap_y(j+2), mx, my)]/2.4e1),
        hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/0.6e1 - mu[idx(i, j, mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.2e1 - mu[idx(i, wrap_y(j+2), mx, my)]/0.6e1),
        hy2inv*(mu[idx(i, j, mx, my)]/0.8e1 - mu[idx(i, wrap_y(j+1), mx, my)]/0.6e1 + mu[idx(i, wrap_y(j+2), mx, my)]/0.8e1)
    ]

    interior_offset = Omega_radius ÷ 2

    # elseif SCHEMA == :mixed_upwind_lambda
    narrow_Omega_radius = narrow_Omega_radius
    narrow_interior_offset = narrow_Omega_radius ÷ 2
    # -- D2, mixed i.e. wide for lambda and narrow for mu --
    for k=1:Omega_radius
        wrapped_x = wrap_x(i - interior_offset - 1 + k)
        wrapped_y = wrap_y(j - interior_offset - 1 + k)
        # -- D1 --
        D1x_centered = 0.5 * (Dpx[k] + Dmx[k])
        D1y_centered = 0.5 * (Dpy[k] + Dmy[k])

        u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
    end
    for k=1:narrow_Omega_radius
        wrapped_x = wrap_x(i - narrow_interior_offset - 1 + k)
        wrapped_y = wrap_y(j - narrow_interior_offset - 1 + k)

        # -- D2, narrow --
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(wrapped_x, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(wrapped_x, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, wrapped_y, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, wrapped_y, mx, my)]
    end

    for kx=1:Omega_radius
        wrapped_kx = wrap_x(i - interior_offset - 1 + kx)

        # Compute inner derivative: ∂x u at position wrapped_kx
        inner_u1_x = 0.0
        inner_u2_x = 0.0
        @inbounds for m=1:Omega_radius
            wrapped_m = wrap_x(wrapped_kx - interior_offset - 1 + m)
            inner_u1_x += Dmx[m] * u1[idx(wrapped_m, j, mx, my)]
            inner_u2_x += Dmx[m] * u2[idx(wrapped_m, j, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = Dpx[kx]
        lambdakx = lambda[idx(wrapped_kx, j, mx, my)]
        mukx = mu[idx(wrapped_kx, j, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += w * lambdakx * inner_u1_x
        # u1_xmux[idx(i, j, mx, my)] += w * mukx * inner_u1_x
        # u2_xmux[idx(i, j, mx, my)] += w * mukx * inner_u2_x
    end

    # y-direction: u1_ymuy = ∂y(mu ∂y u1), u2_ymuy = ∂y(mu ∂y u2), u2_ylambday = ∂y(lambda ∂y u2)
    for ky=1:Omega_radius
        wrapped_ky = wrap_y(j - interior_offset - 1 + ky)

        # Compute inner derivative: ∂y u at position wrapped_ky
        inner_u1_y = 0.0
        inner_u2_y = 0.0
        @inbounds for n=1:Omega_radius
            wrapped_n = wrap_y(wrapped_ky - interior_offset - 1 + n)
            inner_u1_y += Dmy[n] * u1[idx(i, wrapped_n, mx, my)]
            inner_u2_y += Dmy[n] * u2[idx(i, wrapped_n, mx, my)]
        end

        # Apply outer derivative with material coefficients
        w = Dpy[ky]
        lambdaky = lambda[idx(i, wrapped_ky, mx, my)]
        muky = mu[idx(i, wrapped_ky, mx, my)]

        # u1_ymuy[idx(i, j, mx, my)] += w * muky * inner_u1_y
        # u2_ymuy[idx(i, j, mx, my)] += w * muky * inner_u2_y
        u2_ylambday[idx(i, j, mx, my)] += w * lambdaky * inner_u2_y
    end

    for kx=1:Omega_radius
        for ky=1:Omega_radius
            wrapped_x = wrap_x(i - interior_offset - 1 + kx)
            wrapped_y = wrap_y(j - interior_offset - 1 + ky)

            # Centered stencils from averaging upwind operators
            D1x_centered = 0.5 * (Dpx[kx] + Dmx[kx])
            D1y_centered = 0.5 * (Dpy[ky] + Dmy[ky])

            u2_ylambdax[idx(i, j, mx, my)] += Dpx[kx] * lambda[idx(wrapped_x, j, mx, my)] * Dmy[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y_centered * mu[idx(i, wrapped_y, mx, my)] * D1x_centered * u2[idx(wrapped_x, wrapped_y, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += Dpy[ky] * lambda[idx(i, wrapped_y, mx, my)] * Dmx[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x_centered * mu[idx(wrapped_x, j, mx, my)] * D1y_centered * u1[idx(wrapped_x, wrapped_y, mx, my)]
        end
    end

    return nothing
end

function apply_D1_boundary_kernel!(u1_x, u1_y, u2_x, u2_y, u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux, u1, u2, mx, my, 
                            lambda, mu, boundary_radius, Omega_radius, BW_D1,
                            D1x_boundary_W, D1x_boundary_E, D1y_boundary_N, D1y_boundary_S, D1x, D1y)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # 
    # NOTE: Should probably do interior skip here
    #
    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Omega_radius/2)
    interior_offset = Omega_radius ÷ 2
    interior_radius_y = my - 2 * BW_D1
    interior_radius_x = mx - 2 * BW_D1
    if i <= BW_D1
        if j <= BW_D1
            # NW corner!
            # Dx, Dxx, Dy and Dyy at boundary
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u2[idx(k, j, mx, my)]

                u1_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
            for kx=1:boundary_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * lambda[idx(kx, j, mx, my)] * D1y_boundary_N[j, ky] * u2[idx(kx, ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * mu[idx(i, ky, mx, my)] * D1x_boundary_W[i, kx] * u2[idx(kx, ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * lambda[idx(i, ky, mx, my)] * D1x_boundary_W[i, kx] * u1[idx(kx, ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * mu[idx(kx, j, mx, my)] * D1y_boundary_N[j, ky] * u1[idx(kx, ky, mx, my)]
                end
            end
        elseif j >= my - BW_D1 + 1
            # SW corner!
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u2[idx(k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
            for kx=1:boundary_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * lambda[idx(kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u2[idx(kx, my - boundary_radius + ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * mu[idx(i, my - boundary_radius + ky, mx, my)] * D1x_boundary_W[i, kx] * u2[idx(kx, my - boundary_radius + ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * lambda[idx(i, my - boundary_radius + ky, mx, my)] * D1x_boundary_W[i, kx] * u1[idx(kx, my - boundary_radius + ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * mu[idx(kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u1[idx(kx, my - boundary_radius + ky, mx, my)]
                end
            end
        else
            # W edge!
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_W[i, k] * u2[idx(k, j, mx, my)]
            end
            # interior y-derivatives
            for k=1:Omega_radius
                u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                # println("i = $i, j = $j, k = $k, idx(i, j - interior_offset - 1 + k, mx, my) = $(idx(i, j - interior_offset - 1 + k, mx, my))")
            end
            for kx=1:boundary_radius
                for ky=1:Omega_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * lambda[idx(kx, j, mx, my)] * D1y[ky] * u2[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_boundary_W[i, kx] * u2[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_boundary_W[i, kx] * u1[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_W[i, kx] * mu[idx(kx, j, mx, my)] * D1y[ky] * u1[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                end
            end
        end
    elseif i >= mx - BW_D1 + 1
        if j <= BW_D1
            # NE corner
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]

                u1_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
            for kx=1:boundary_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * lambda[idx(mx - boundary_radius + kx, j, mx, my)] * D1y_boundary_N[j, ky] * u2[idx(mx - boundary_radius + kx, ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * mu[idx(i, ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u2[idx(mx - boundary_radius + kx, ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * lambda[idx(i, ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u1[idx(mx - boundary_radius + kx, ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * mu[idx(mx - boundary_radius + kx, j, mx, my)] * D1y_boundary_N[j, ky] * u1[idx(mx - boundary_radius + kx, ky, mx, my)]
                end
            end
        elseif j >= my - BW_D1 + 1
            # SE corner
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
            for kx=1:boundary_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * lambda[idx(mx - boundary_radius + kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u2[idx(mx - boundary_radius + kx, my - boundary_radius + ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * mu[idx(i, my - boundary_radius + ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u2[idx(mx - boundary_radius + kx, my - boundary_radius + ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * lambda[idx(i, my - boundary_radius + ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u1[idx(mx - boundary_radius + kx, my - boundary_radius + ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * mu[idx(mx - boundary_radius + kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u1[idx(mx - boundary_radius + kx, my - boundary_radius + ky, mx, my)]
                end
            end
        else
            # E edge
            for k=1:boundary_radius
                u1_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]
            end
            # interior y-derivatives
            for k=1:Omega_radius
                u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
            end
            for kx=1:boundary_radius
                for ky=1:Omega_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * lambda[idx(mx - boundary_radius + kx, j, mx, my)] * D1y[ky] * u2[idx(mx - boundary_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u2[idx(mx - boundary_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * u1[idx(mx - boundary_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x_boundary_E[i - (BW_D1 + interior_radius_x), kx] * mu[idx(mx - boundary_radius + kx, j, mx, my)] * D1y[ky] * u1[idx(mx - boundary_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                end
            end
        end
    else
        if j <= BW_D1
            # N edge
            for k=1:boundary_radius
                u1_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
            # interior x-derivatives
            for k=1:Omega_radius
                u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
            # Mixed derivatives
            for kx=1:Omega_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x[kx] * lambda[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_boundary_N[j, ky] * u2[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * mu[idx(i, ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_N[j, ky] * lambda[idx(i, ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x[kx] * mu[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_boundary_N[j, ky] * u1[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                end
            end
        elseif j >= my - BW_D1 + 1
            # S edge
            for k=1:boundary_radius
                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
            # interior x-derivatives
            for k=1:Omega_radius
                u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
            # Mixed derivatives
            for kx=1:Omega_radius
                for ky=1:boundary_radius
                    u2_ylambdax[idx(i, j, mx, my)] += D1x[kx] * lambda[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u2[idx(i - interior_offset - 1 + kx, my - boundary_radius + ky, mx, my)]
                    u2_xmuy[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * mu[idx(i, my - boundary_radius + ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, my - boundary_radius + ky, mx, my)]
                    u1_xlambday[idx(i, j, mx, my)] += D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * lambda[idx(i, my - boundary_radius + ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, my - boundary_radius + ky, mx, my)]
                    u1_ymux[idx(i, j, mx, my)] += D1x[kx] * mu[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_boundary_S[j - (BW_D1 + interior_radius_y), ky] * u1[idx(i - interior_offset - 1 + kx, my - boundary_radius + ky, mx, my)]
                end
            end
        end
    end

    return nothing
end

function apply_D1_inner_kernel!(u1_x, u1_y, u2_x, u2_y, u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                                u1, u2, mx, my, Omega_radius, lambda, mu, D1x, D1y, BW_D1)
    
    i = BW_D1 + (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = BW_D1 + (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i >= mx - BW_D1 + 1) || (j >= my - BW_D1 + 1) || (i <= BW_D1) || (j <= BW_D1)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    interior_offset = Omega_radius ÷ 2
    # Interior points
    for k=1:Omega_radius
        u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
    end
    for kx=1:Omega_radius
        for ky=1:Omega_radius
            u2_ylambdax[idx(i, j, mx, my)] += D1x[kx] * lambda[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y[ky] * u2[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x[kx] * mu[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y[ky] * u1[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
        end
    end

    return nothing
end

function apply_D2_narrow_boundary_kernel!(u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux, u1, u2, 
    mx, my, lambda, mu, boundary_radius, Omega_radius, BW_D2, hx, hy,
    D2x_lambda_boundary_W_vec, D2x_mu_boundary_W_vec, D2x_lambda_boundary_E_vec, D2x_mu_boundary_E_vec,
    D2y_lambda_boundary_N_vec, D2y_mu_boundary_N_vec, D2y_lambda_boundary_S_vec, D2y_mu_boundary_S_vec,
    D2x_lambda_vec, D2x_mu_vec, D2y_lambda_vec, D2y_mu_vec)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i > mx) || (j > my)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    # -- NOTE: Kept for now since it might be faster when AI > BW --
    # hx2inv = 1.0f0/(hx*hx)
    # hy2inv = 1.0f0/(hy*hy)
    # @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    # @inline wrap_y(idx) = (idx+my-1) % my + 1
 
    # D2x_lambda = -@SVector [
    #     hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
    #     hx2inv*(-lambda[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
    #     hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + lambda[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
    #     hx2inv*(-lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
    #     hx2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    # ]
    # D2x_mu = -@SVector [
    #     hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
    #     hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
    #     hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + mu[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
    #     hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
    #     hx2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + mu[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    # ]
    # D2y_mu = -@SVector [
    #     hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
    #     hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
    #     hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + mu[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
    #     hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
    #     hy2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + mu[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    # ]
    # D2y_lambda = -@SVector [
    #     hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
    #     hy2inv*(-lambda[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
    #     hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + lambda[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
    #     hy2inv*(-lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
    #     hy2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    # ]

    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Omega_radius/2)
    interior_offset = Omega_radius ÷ 2
    interior_radius_y = my - 2 * BW_D2
    interior_radius_x = mx - 2 * BW_D2

    if i <= BW_D2
        # Load West boundary stencils
        D2x_lambda_boundary_W = D2x_lambda_boundary_W_vec[j]
        D2x_mu_boundary_W = D2x_mu_boundary_W_vec[j]
        
        if j <= BW_D2
            # NW corner - also need North boundary stencils
            D2y_lambda_boundary_N = D2y_lambda_boundary_N_vec[i]
            D2y_mu_boundary_N = D2y_mu_boundary_N_vec[i]
            # Dx, Dxx, Dy and Dyy at boundary
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u2[idx(k, j, mx, my)]

                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u2[idx(i, k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
        elseif j >= my - BW_D2 + 1
            # SW corner - also need South boundary stencils
            D2y_lambda_boundary_S = D2y_lambda_boundary_S_vec[i]
            D2y_mu_boundary_S = D2y_mu_boundary_S_vec[i]
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u2[idx(k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
        else
            # W edge - need interior y stencils
            D2y_mu = D2y_mu_vec[i, j]
            D2y_lambda = D2y_lambda_vec[i, j]
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_W[i, k] * u2[idx(k, j, mx, my)]
            end
            # interior y-derivatives
            for k=1:Omega_radius
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                # println("i = $i, j = $j, k = $k, idx(i, j - interior_offset - 1 + k, mx, my) = $(idx(i, j - interior_offset - 1 + k, mx, my))")
            end
        end
    elseif i >= mx - BW_D2 + 1
        # Load East boundary stencils
        D2x_lambda_boundary_E = D2x_lambda_boundary_E_vec[j]
        D2x_mu_boundary_E = D2x_mu_boundary_E_vec[j]
        
        if j <= BW_D2
            # NE corner - also need North boundary stencils
            D2y_lambda_boundary_N = D2y_lambda_boundary_N_vec[i]
            D2y_mu_boundary_N = D2y_mu_boundary_N_vec[i]
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]

                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u2[idx(i, k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
        elseif j >= my - BW_D2 + 1
            # SE corner - also need South boundary stencils
            D2y_lambda_boundary_S = D2y_lambda_boundary_S_vec[i]
            D2y_mu_boundary_S = D2y_mu_boundary_S_vec[i]
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
        else
            # E edge - need interior y stencils
            D2y_mu = D2y_mu_vec[i, j]
            D2y_lambda = D2y_lambda_vec[i, j]
            for k=1:boundary_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u1[idx(mx - boundary_radius + k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu_boundary_E[i - (BW_D2 + interior_radius_x), k] * u2[idx(mx - boundary_radius + k, j, mx, my)]
            end
            # interior y-derivatives
            for k=1:Omega_radius
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
            end
        end
    else
        if j <= BW_D2
            # N edge - need North boundary and interior x stencils
            D2y_lambda_boundary_N = D2y_lambda_boundary_N_vec[i]
            D2y_mu_boundary_N = D2y_mu_boundary_N_vec[i]
            D2x_lambda = D2x_lambda_vec[i, j]
            D2x_mu = D2x_mu_vec[i, j]
            for k=1:boundary_radius
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u1[idx(i, k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_N[j, k] * u2[idx(i, k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_N[j, k] * u2[idx(i, k, mx, my)]
            end
            # interior x-derivatives
            for k=1:Omega_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
        elseif j >= my - BW_D2 + 1
            # S edge - need South boundary and interior x stencils
            D2y_lambda_boundary_S = D2y_lambda_boundary_S_vec[i]
            D2y_mu_boundary_S = D2y_mu_boundary_S_vec[i]
            D2x_lambda = D2x_lambda_vec[i, j]
            D2x_mu = D2x_mu_vec[i, j]
            for k=1:boundary_radius
                # Since stencils already are mirrored, we go from down to up
                u1_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u1[idx(i, my - boundary_radius + k, mx, my)]
                u2_ymuy[idx(i, j, mx, my)] += D2y_mu_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
                u2_ylambday[idx(i, j, mx, my)] += D2y_lambda_boundary_S[j - (BW_D2 + interior_radius_y), k] * u2[idx(i, my - boundary_radius + k, mx, my)]
            end
            # interior x-derivatives
            for k=1:Omega_radius
                u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
        end
    end

    return nothing
end

function apply_D2_narrow_inner_kernel!(u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux, u1, u2, 
    mx, my, lambda, mu, Omega_radius, BW_D2, hx, hy, D2x_lambda_vec, D2x_mu_vec, D2y_lambda_vec, D2y_mu_vec)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i >= mx - BW_D2 + 1) || (j >= my - BW_D2 + 1) || (i <= BW_D2) || (j <= BW_D2)
        return nothing
    end

    D2x_lambda = D2x_lambda_vec[i, j]
    D2x_mu = D2x_mu_vec[i, j]

    D2y_mu = D2y_mu_vec[i, j]
    D2y_lambda = D2y_lambda_vec[i, j]

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    # -- NOTE: Kept for now since it might be faster when AI > BW --
    # # TODO: Remove wrap from all kernels
    # @inline wrap_x(idx) = (idx+mx-1) % mx + 1
    # @inline wrap_y(idx) = (idx+my-1) % my + 1
    
    # hx2inv = 1.0f0/(hx*hx)
    # D2x_lambda = -@SVector [
    #     hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
    #     hx2inv*(-lambda[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - lambda[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
    #     hx2inv*(lambda[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + lambda[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
    #     hx2inv*(-lambda[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - lambda[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
    #     hx2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + lambda[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    # ]
    # D2x_mu = -@SVector [
    #     hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/8.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
    #     hx2inv*(-mu[idx(wrap_x(i-2), j, mx, my)]/6.0f0 - mu[idx(wrap_x(i-1), j, mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0),
    #     hx2inv*(mu[idx(wrap_x(i-2), j, mx, my)]/24.0f0 + mu[idx(wrap_x(i-1), j, mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(wrap_x(i+1), j, mx, my)]/1.2f0 + mu[idx(wrap_x(i+2), j, mx, my)]/24.0f0),
    #     hx2inv*(-mu[idx(wrap_x(i-1), j, mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/2.0f0 - mu[idx(wrap_x(i+2), j, mx, my)]/6.0f0),
    #     hx2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(wrap_x(i+1), j, mx, my)]/6.0f0 + mu[idx(wrap_x(i+2), j, mx, my)]/8.0f0)
    # ]
    # hy2inv = 1.0f0/(hy*hy)
    # D2y_mu = -@SVector [
    #     hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + mu[idx(i, j, mx, my)]/8.0f0),
    #     hy2inv*(-mu[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - mu[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
    #     hy2inv*(mu[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + mu[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + mu[idx(i, j, mx, my)]*0.75f0 + mu[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + mu[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
    #     hy2inv*(-mu[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - mu[idx(i, j, mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - mu[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
    #     hy2inv*(mu[idx(i, j, mx, my)]/8.0f0 - mu[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + mu[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    # ]
    # D2y_lambda = -@SVector [
    #     hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 + lambda[idx(i, j, mx, my)]/8.0f0),
    #     hy2inv*(-lambda[idx(i, wrap_y(j-2), mx, my)]/6.0f0 - lambda[idx(i, wrap_y(j-1), mx, my)]/2.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0),
    #     hy2inv*(lambda[idx(i, wrap_y(j-2), mx, my)]/24.0f0 + lambda[idx(i, wrap_y(j-1), mx, my)]/1.2f0 + lambda[idx(i, j, mx, my)]*0.75f0 + lambda[idx(i, wrap_y(j+1), mx, my)]/1.2f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/24.0f0),
    #     hy2inv*(-lambda[idx(i, wrap_y(j-1), mx, my)]/6.0f0 - lambda[idx(i, j, mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/2.0f0 - lambda[idx(i, wrap_y(j+2), mx, my)]/6.0f0),
    #     hy2inv*(lambda[idx(i, j, mx, my)]/8.0f0 - lambda[idx(i, wrap_y(j+1), mx, my)]/6.0f0 + lambda[idx(i, wrap_y(j+2), mx, my)]/8.0f0)
    # ]

    interior_offset = Omega_radius ÷ 2

    for k=1:Omega_radius
        u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_ylambday[idx(i, j, mx, my)] += D2y_lambda[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
    end

    return nothing
end

# Fuse apply_D1_inner_kernel and apply_D2_narrow_inner_kernel
function apply_narrow_inner_kernel!(u1_x, u1_y, u2_x, u2_y, u2_ylambdax, u2_xmuy, u1_xlambday, u1_ymux,
                                u1, u2, mx, my, Omega_radius, lambda, mu, D1x, D1y, BW_D1,
                                u1_xlambdax, u1_xmux, u1_ymuy, u2_ylambday, u2_ymuy, u2_xmux,
                                BW_D2, hx, hy, D2x_lambda_vec, D2x_mu_vec, D2y_lambda_vec, D2y_mu_vec)
    i = BW_D1 + (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = BW_D1 + (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i >= mx - BW_D1 + 1) || (j >= my - BW_D1 + 1) || (i <= BW_D1) || (j <= BW_D1)
        return nothing
    end

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    D2x_lambda = D2x_lambda_vec[i, j]
    D2x_mu = D2x_mu_vec[i, j]

    D2y_mu = D2y_mu_vec[i, j]
    D2y_lambda = D2y_lambda_vec[i, j]

    interior_offset = Omega_radius ÷ 2
    # Interior points
    for k=1:Omega_radius
        u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]

        u1_xlambdax[idx(i, j, mx, my)] += D2x_lambda[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u1_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_xmux[idx(i, j, mx, my)] += D2x_mu[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_ymuy[idx(i, j, mx, my)] += D2y_mu[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_ylambday[idx(i, j, mx, my)] += D2y_lambda[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
        
        for ky=1:Omega_radius
            u2_ylambdax[idx(i, j, mx, my)] += D1x[k] * lambda[idx(i - interior_offset - 1 + k, j, mx, my)] * D1y[ky] * u2[idx(i - interior_offset - 1 + k, j - interior_offset - 1 + ky, mx, my)]
            u2_xmuy[idx(i, j, mx, my)] += D1y[ky] * mu[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[k] * u2[idx(i - interior_offset - 1 + k, j - interior_offset - 1 + ky, mx, my)]
            u1_xlambday[idx(i, j, mx, my)] += D1y[ky] * lambda[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[k] * u1[idx(i - interior_offset - 1 + k, j - interior_offset - 1 + ky, mx, my)]
            u1_ymux[idx(i, j, mx, my)] += D1x[k] * mu[idx(i - interior_offset - 1 + k, j, mx, my)] * D1y[ky] * u1[idx(i - interior_offset - 1 + k, j - interior_offset - 1 + ky, mx, my)]
        end
    end

    return nothing
end