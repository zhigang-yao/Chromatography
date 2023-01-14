function [L,G]=Grad(y_obs,Tn,S,R)
% Function to calculate the numerical gradient
    h = 0.000001;
    Ss=[S;S+[h,0];S-[h,0];S+[0,h];S-[0,h]];
    xi = recover(R,Ss);
    Y = solver(xi,Tn);
    E = Y - y_obs;
    Ls = diag(E*E')';
    L = Ls(1);
    G = [Ls(2)-Ls(3),Ls(4)-Ls(5)]/(2*h);
end