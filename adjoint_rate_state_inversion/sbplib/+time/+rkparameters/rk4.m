function [a,b,c,s] = rk4()

% Butcher tableau for classical RK$
s = 4;
a = sparse(s,s);
a(2,1) = 1/2;
a(3,2) = 1/2;
a(4,3) = 1;
b = 1/6*[1; 2; 2; 1];
c = [0; 1/2; 1/2; 1];

end