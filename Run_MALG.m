%% Run MALG for N times
clear
load data.mat

K = 10000;
N = 100;
AR = zeros(N,3);
[sigma2,r1,r2,s1,s2,xi1,xi2,xi3,xi4] = deal(zeros(N,K));

parfor itera = 1:N
    BN=500;
    [xi,sig2,Rs,S,ar]=MALG(y_obs, Tn, 10, K+BN,'Loss-Ratio.mat');
    
    xi = xi(BN+1:end,:);
    Rs = Rs(BN+1:end,:);
    S = S(BN+1:end,:);
    sig2 = sig2(BN+1:end,:);
    
    AR(itera,:) = ar';
    r1(itera,:) = Rs(:,1)';
    r2(itera,:) = Rs(:,2)';
    s1(itera,:) = S(:,1)';
    s2(itera,:) = S(:,2)';
    sigma2(itera,:) = sig2;
    xi1(itera,:) = xi(:,1)';
    xi2(itera,:) = xi(:,2)';
    xi3(itera,:) = xi(:,3)';
    xi4(itera,:) = xi(:,4)';
end
save MALG_GMM_2C.mat AR sigma2 r1 r2 s1 s2 xi1 xi2 xi3 xi4

