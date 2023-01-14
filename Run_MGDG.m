%% Run MGDG for N times
clear
load data.mat

K = 10000;
N = 100;
AR = zeros(N,2);

[sigma2,r1,r2,xi1,xi2,xi3,xi4,xi_r1,xi_r2,xi_r3,xi_r4] = deal(zeros(N,K));


parfor itera = 1:N
    BN=500;
    [xi,sig2,Rs,ar]=MGDG(y_obs, Tn, 10, K+BN,'Loss-Ratio.mat');
    
    xi = xi(BN+1:end,:);
    Rs = Rs(BN+1:end,:);
    sig2 = sig2(BN+1:end,:);
    xi_rec = xi;
    
    for i = 2:K
        Rm = mean(Rs(1:i,:));
        [S,L] = GD2D(y_obs,Tn,[3,3],Rm);
        xi_r = recover(Rm,S);
        xi_rec(i,:) = xi_r;
    end
    
    AR(itera,:) = ar';
    r1(itera,:) = Rs(:,1)';
    r2(itera,:) = Rs(:,2)';
    sigma2(itera,:) = sig2;
    xi1(itera,:) = xi(:,1)';
    xi2(itera,:) = xi(:,2)';
    xi3(itera,:) = xi(:,3)';
    xi4(itera,:) = xi(:,4)';
    
    xi_r1(itera,:) = xi_r(:,1)';
    xi_r2(itera,:) = xi_r(:,2)';
    xi_r3(itera,:) = xi_r(:,3)';
    xi_r4(itera,:) = xi_r(:,4)';

end
save MGDG_GMM_2C.mat sigma2 r1 r2 xi1 xi2 xi3 xi4 xi_r1 xi_r2 xi_r3 xi_r4 AR


