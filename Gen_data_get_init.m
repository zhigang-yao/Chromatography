% Script for generate date and collect information for initilization

addpath(genpath(pwd));
%% Generate Data
clear
clc

xi_t = recover([1/3,2/3],[1,4]);
fprintf("ground truth: (%.4f,%.4f,%.4f,%.4f)\n",xi_t);
sigma2_t = 0.001;
Tn = linspace(-2,7,200);

Y  = solver(xi_t,Tn);
noise = sqrt(sigma2_t) * randn(1,size(Tn,2));
y_obs = Y + noise;

plot(Tn,[y_obs;Y]);

save data.mat y_obs Tn xi_t

%% Find init
Size = 1000;

Ras = gallery('uniformdata',[1,Size],0);
Rbs = gallery('uniformdata',[1,Size],1000);

mSas = zeros(1,Size);
mSbs = zeros(1,Size);
mLs  = zeros(1,Size);


S0 = [2,2];
parfor i = 1:Size
    R = [Ras(i),Rbs(i)];
    [S,L] = GD2D(y_obs,Tn,S0,R)
    mSas(i) = S(1);
    mSbs(i) = S(2);
    mLs(i)  = L;
end
[xq,yq] = meshgrid(linspace(0.02,0.98,Size/2), linspace(0.02,0.98,Size/2));
v1q = griddata(Ras,Rbs,mLs,xq,yq);
mesh(xq,yq,v1q)
hold on
plot3(Ras,Rbs,mLs,'o','Color',[0.8500 0.3250 0.0980]','MarkerSize',2)
xlabel('\eta_1')
ylabel('\eta_2')
zlabel('Loss')
view(-30,11)

ax1 = axes('Position',[0.7 0.65 0.25 0.25]);
axes(ax1)
mesh(xq,yq,v1q)
xlim([0,1])
ylim([0,1])
view(0,90)
hold off

save Loss-ratio.mat Ras Rbs mLs mSas mSbs Size



