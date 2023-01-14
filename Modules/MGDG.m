function [xi,sigma2,Rs,AccRate]=MGDG(Y_obs, Tn, phi, K,InputName)
% Func: preform MGDG sampling
    %  Input: Y_obs - m x N matrix, observed data
    %         Tn    - 1 x N, corresponding time
    %         phi   - Some hyper-parameter
    %         K     - Number of samples required to be sampled
    %         m     - Number of observed y_obs
    %         n     - number of linspace to be used
    %     InputName - InputName.mat file will be loaded for initialization  
    %  Output: xi   - Parameter of interest
    %       sigma2  - Noise variance
    %         Rs    - Compressed parameter R
    %     AccRate   - Acceptance rate from each parameter
  %% Calculate Prior
  load(InputName, 'Ras', 'Rbs', 'mLs')
  [xq,yq] = meshgrid(linspace(0,1,length(Ras)/2), linspace(0,1,length(Ras)/2));
  vq = exp(- 8 * mLs);
  vq = griddata(Ras,Rbs,vq,xq,yq);
  AccCount = zeros(2,2);

  %% Initialization   
    dem = 2;
    tn = Tn;
    n = size(tn,2);
    m=length(Y_obs(:,1));    
    y_obs=zeros(m,n);    
    for j=1:m        
        y_obs(j,:)=interp1(Tn,Y_obs(j,:),tn);
    end
    
    xi= zeros(K,4);
    sigma2= zeros(K,1);
    Rs = zeros(K,2);
    
    Rs_old = [Ras(mLs == min(mLs)),Rbs(mLs == min(mLs))];
    [xi_old,S_old] = Get_xi(Y_obs,[2,2],Rs_old,Tn);
    phi_sig = (min(mLs)/size(Tn,2));
    
    a_old = rpi_a(phi);
    sigma2_old = rpi_sig(phi_sig);
    sig_a0 = 1;
    sig_b0 = 3*phi_sig;

    y_old = solver(xi_old,tn);    
    e_old = y_old-y_obs;
    eTe_old = e_old'*e_old;

    M_old=Kernel(a_old,tn)+sigma2_old*eye(n);
    M_old_inv=inv(M_old);
    M_old_det=det(M_old);
 
  %% Sampling iterations
    for i=1:(K)
      %% Update sigam2
        sig_a = sig_a0 + n/2;
        sig_b = sig_b0 + trace(eTe_old)/2;
        sigma2_old = 1 / (gamrnd(sig_a , 1/sig_b));
        M_new = sigma2_old*eye(n);
        M_old_inv=inv(M_new);
        M_old_det=det(M_new);
        
      %% Update a
          a_new = a_old;
          a_old = a_new;

      %% Update R
        Rs_new = Rs_old;
        for j = 1:dem
            AccCount(j,2) = AccCount(j,2) + 1;
            Cand = Rs_old;
            Rs_0 = Rs_old(j);
            Rs_1 = rq_R(Rs_0);
            
            Cand(j) = Rs_1;
            [xi_cand,S_cand] = Get_xi(Y_obs,[2,2],Cand,Tn);
            y_1 = solver(xi_cand,tn);
            e_1 = y_1-y_obs;
            eTe_1 = e_1'*e_1;
            
            eTe_0 = eTe_old;
            M_inv = M_old_inv;
            
            logrho1 = trace((eTe_0-eTe_1)*M_inv)/2;
            logrho2 = log(interp2(xq,yq,vq,Cand(1),Cand(2)))-log(interp2(xq,yq,vq,Rs_old(1),Rs_old(2)));
            logrho3 = log(dq_R(Rs_0,Rs_1))-log(dq_R(Rs_1,Rs_0));
            logrho  = logrho1+logrho2+logrho3;
            
            alpha01 = exp(min(0,logrho));
            if log(rand())<min(0,logrho)
                Rs_new(j) = Rs_1;
                AccCount(j,1) = AccCount(j,1) + 1;
            end
        end
        Rs_old  = Rs_new;
        [xi_old,~] = Get_xi(Y_obs,[2,2],sort(Rs_old),Tn);
        y_old = solver(xi_old,tn);    
        e_old = y_old-y_obs;
        eTe_old = e_old'*e_old;
        xi(i,:)=xi_old;
        sigma2(i,:)=sigma2_old;
        Rs(i,:) = sort(Rs_old);
        AccRate = AccCount(:,1)./AccCount(:,2);
        
    end
end


function k=Kernel(a,Tn)
    k=exp(-(Tn-Tn').^2/(2*a^2))*0;
end


function tn = rtnorm(mean,sd)
    tn = normrnd(mean,sd);
    while tn<=0.00001
        tn = normrnd(mean,sd);
    end
end
function qtn = dtnorm(x,mean,sd)
    qtn = normpdf(x,mean,sd)/(1-normcdf(0,mean,sd));
end

%% Functions: pi()
function [a] = rpi_a(phi)
    a = rtnorm(phi,4);
end

function [sigma2] = rpi_sig(phi_sig)
    sigma2 = 1/gamrnd(1,1/(3*phi_sig));
end

%% Functions: q()
function [R_new] = rq_R(R_old)
    R_new = rtnorm(R_old,0.02);
    R_new = R_new - 2*(R_new - 0.99)*(R_new>0.99);
    R_new = R_new + 2*(0.01 - R_new) *(R_new < 0.01);
end
function [R_new] = dq_R(R_new,R_old)
    R_new = dtnorm(R_new,R_old,0.02);
end

function [xi,S] = Get_xi(y_obs,S,R,Tn)
    [S,~] = GD2D(y_obs,Tn,S,R);
    xi = recover(R,S);
end