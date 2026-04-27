function A = cell_row_to_diag_blockmatrix(A_cell)
    n = numel(A_cell);
    A = cell(n,n);
    for i = 1:n
        A{i,i} = A_cell{i};
    end
    A = blockmatrix.toMatrix(A);
end