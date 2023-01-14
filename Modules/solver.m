function y = solver(xi,Tn)
% Solver handling xi with matrix calculation
     s = [xi(:,1)+xi(:,2),xi(:,3)+xi(:,4)];
     r = [xi(:,1),xi(:,3)]./ s ;
     
%% Par solver
    y = r(:,1).* exp(-(Tn - s(:,1)).^2/(2))/(2*pi)^.5 + r(:,2) .*exp(-(Tn - s(:,2)).^2/(2))/(2*pi)^.5;
    
end

