function y = solver_xi(xi,Tn)
% solver with xi
r = compress(xi);
s = [xi(1)+xi(2),xi(3)+xi(4)];
y = r(1) * exp(-(Tn - s(1)).^2/(2))/(2*pi)^.5 + r(2) *exp(-(Tn - s(2)).^2/(2))/(2*pi)^.5;
    
end

