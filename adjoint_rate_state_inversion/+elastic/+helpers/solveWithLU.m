function x = solveWithLU(A, b)

	x = A.Q*(A.U\(A.L\(A.P*b)));

end