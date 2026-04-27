function inds = d2_sparsity_pattern_inds(m, order, BP, interior_offset, closure_offset)
    sparsity_pattern = sparse(m,m);
    inner_stencil_inds = -order/2-interior_offset:order/2+interior_offset;
    inner_stencil_pattern = ones(m,length(inner_stencil_inds));
    sparsity_pattern = spdiags(inner_stencil_pattern,inner_stencil_inds,sparsity_pattern);
    sparsity_pattern(1:BP,1:BP+closure_offset) = 1;
    sparsity_pattern(end-BP+1:end,end-(BP+closure_offset)+1:end) = 1;
    for k = 1:closure_offset
        sparsity_pattern(BP+k,1:BP+closure_offset+k) = 1;
        sparsity_pattern(end-(BP+k)+1,end-(BP+k+closure_offset)+1:end) = 1;
    end
    inds  = find(sparsity_pattern);
end