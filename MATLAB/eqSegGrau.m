function [x1,x2] = eqSegGrau(a,b,c)
    delta = b^2-4*a*c;
    sqrt_delta = sqrt(delta);
    x1 = (-b + sqrt(delta))/(2*a);
    x2 = (-b - sqrt(delta))/(2*a);
end