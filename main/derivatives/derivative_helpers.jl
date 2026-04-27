@inline function extract_row_x!(buffer, material_field, j, mx, my)
    @inbounds for i in 1:mx
        buffer[i] = material_field[idx(i, j, mx, my)]
    end
    return buffer
end

@inline function extract_col_y!(buffer, material_field, i, mx, my)
    @inbounds for j in 1:my
        buffer[j] = material_field[idx(i, j, mx, my)]
    end
    return buffer
end