# idx for (i,j) in 2D grid stored in 1D array
idx(i, j, mx, my) = (j-1)*mx + i

# # Indexing functions for staggered grid variables with different dimensions
# # v1 at (primal_x, dual_y): mx × (my+1) for non-periodic
# idx_v1(i, j, mx, my, periodic) = idx(i, j, mx, my)
# # v2 at (dual_x, primal_y): (mx+1) × my for non-periodic
# idx_v2(i, j, mx, my, periodic) = periodic ? idx(i, j, mx, my) : (j-1)*(mx+1) + i
# # σ11, σ22 at (dual_x, dual_y): (mx+1) × (my+1) for non-periodic
# idx_sigma_dual(i, j, mx, my, periodic) = periodic ? idx(i, j, mx, my) : (j-1)*(mx+1) + i
# # σ12 at (primal_x, primal_y): mx × my
# idx_sigma12(i, j, mx, my, periodic) = idx(i, j, mx, my)