function [xi,sigma2]=MMH(Y_obs, Tn, phi, K)
    %  Input: Y_obs - m x N matrix, observed data
    %     Injections- m x 2 matrix, preseted parameters
    %         Tn    - 1 x N, corresponding time
    %         phi   - Some hyper-parameter
    %         K     - Number of samples required to be sampled
    %         m     - Number of observed y_obs
    %         n     - number of linspace to be used
    %  Output:xi    - parameter of interest
    %         sigma2- noise variance
    

  %% Initialization   
    dem = 4;
    m=length(Y_obs(:,1));   
    tn = Tn;
    n = size(tn,2); 

    y_obs=zeros(m,n);    
    for j=1:m        
        y_obs(j,:)=interp1(Tn,Y_obs(j,:),tn);
    end
    
    xi= zeros(K,dem);
    sigma2= zeros(K,1);
 
    xi_old = zeros(1,dem);
    for i = 1:dem        
        xi_old(i) = rpi_xi(phi);
    end
    xi_old = [0.5,0.5,3,1];
    
    a_old = rpi_a(phi);
    sig_a0 = 1;
    sig_b0 = 3*0.02;
    sigma2_old = rpi_sig(sig_b0);
    sigma2_old = 0.001;
    phi_sig = 0.001;

    y_old = solver_xi(xi_old,tn);   
    e_old = y_old-y_obs;
    eTe_old=e_old'*e_old;

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

     %% Update xi
        xi_new = xi_old;
        for j = 1:length(xi_old)
            Cand = xi_old;
            cand = rq_xi(xi_old(j));
            Cand(j) = cand;

            y_new = solver_xi(Cand,tn);
            e_new=y_new-y_obs;
            eTe_new=e_new'*e_new;
            
            logrho1 = trace((eTe_old-eTe_new)*M_old_inv)/2;
            logrho2 = log(dpi_xi(cand,phi))-log(dpi_xi(xi_old(j),phi));
            logrho3 = log(dq_xi(xi_old(j),cand))-log(dq_xi(cand,xi_old(j)));
            logrho  = logrho1+logrho2+logrho3;
            
            if log(rand())<min(0,logrho)
                xi_new(j) = cand;
            end
        end

        xi_old = xi_new;
        y_old = solver_xi(xi_old,tn);    
        e_old = y_old-y_obs;
        eTe_old = e_old'*e_old;
        
        xi(i,:)=xi_old;
        sigma2(i,:)=sigma2_old;
    end
end


function k=Kernel(a,Tn)
    k=exp(-(Tn-Tn').^2/(2*a^2))*0;
end

function tn = rtnorm(mean,sd)
    tn = normrnd(mean,sd);
    while tn<=0
        tn = normrnd(mean,sd);
    end
end
function qtn = dtnorm(x,mean,sd)
    qtn = normpdf(x,mean,sd)/(1-normcdf(0,mean,sd));
end

%% Functions: pi() for 3 parameters
function [a] = rpi_a(phi)
    a = rtnorm(phi,4);
end

function [sigma2] = rpi_sig(phi_sig)
    sigma2 = 1/gamrnd(1,1/(3*phi_sig));
end
function [d] = dpi_sig(sigma2,phi_sig)
    d = gampdf(1/sigma2,1,1/(3*phi_sig));
end

function [xi] = rpi_xi(phi)
    xi = rtnorm(phi/10,3);
end
function [d] = dpi_xi(xi,phi)
    d = dtnorm(xi,phi/10,3);
end
%% Functions: q()

function [xi_new] = rq_xi(xi_old)
    xi_new = rtnorm(xi_old,0.05);
end
function [xi_new] = dq_xi(xi_new,xi_old)
    xi_new = dtnorm(xi_new,xi_old,0.05);
end

