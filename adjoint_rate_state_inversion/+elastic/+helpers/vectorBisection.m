function x = vectorBisection(f, interval, tol)

	a = interval(:,1);
	b = interval(:,2);
	x = (a+b)/2;
	err = max(b-a)/2;

	fa = f(a);

	while err > tol

		fx = f(x);

		% Find sign change
		ind = fa.*fx < 0;

		% Set new limits and guess
		a = ind.*a + (~ind).*x;
		b = ind.*x + (~ind).*b;
		x = (a+b)/2;

		% Update left function values
		fa = ind.*fa + (~ind).*fx;

		err = err/2;

	end

end