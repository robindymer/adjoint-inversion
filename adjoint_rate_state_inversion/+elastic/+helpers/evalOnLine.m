function fval = evalOnLine(X, f)
	if isa(f, 'function_handle')
		fval = f(X(:,1), X(:,2));
	else
		% Assume that f is a scalar
		fval = 0*X(:,1) + f;
	end
end