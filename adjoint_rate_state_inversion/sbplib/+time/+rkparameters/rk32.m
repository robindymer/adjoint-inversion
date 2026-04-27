function [a,b,bStar,c,s] = rk32()

% Butcher tableau for embedded adaptive RK3-2
s = 3;
a = sparse(s,s);
a(2,1) = 1/2;
a(3,1) = -1;
a(3,2) = 2;
b = 1/6*[1; 4; 1];
bStar = 1/2*[1; 0; 1];
c = [0; 1/2; 1];

end