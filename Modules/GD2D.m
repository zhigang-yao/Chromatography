function [S,L] = GD2D(y_obs,Tn,S0,R)
% Function: gradient descent in 2D case
    alpha = 0.5;
    S = S0;
    
    for ite = 1:1000   
        [L,G]=Grad(y_obs,Tn,S,R);
        if norm(G) < 0.0001
           break
        end
        Sn = S - alpha * G / norm(G);        
        t = 0;
        while Loss(y_obs,Tn,Sn,R) > L
            t = t+1;
            Sn = S - 0.9^t * alpha * G / norm(G);
            if  0.9^t * alpha < 0.00001
                break
            end
        end
        
        S = Sn;
        S(S<0)=0.002;
    end
    if ite == 1000
        fprintf("Warning from GD2D:,%.4f,%.4f\n",R);
    end
end

function [L]=Loss(y_obs,Tn,S,R)
        xi = recover(R,S);
        y = solver(xi,Tn);
        e = y_obs - y;
        L = e*e';
end
