clear
load('MGDG_GMM_2C.mat');
figure(1)
rs = [r1(1,:);r2(1,:)]';
r_t = [1/3,2/3];
es =  (rs-r_t);
[S1,AX1,BigAx1,H1,HAx1] = plotmatrix(es);
DrawPNG('Simu_case_1_PlotMatrix_r');

figure(2)
xis = [xi1(1,:);xi2(1,:);xi3(1,:);xi4(1,:)]';
xi_t = [1/3,2/3,8/3,4/3];
xis = (xis-xi_t);
[S2,AX2,BigAx2,H2,HAx2] = plotmatrix(xis);
DrawPNG('Simu_case_1_PlotMatrix_xi');

inds = (1:20)*500;

vars = {'r1','r2','xi1','xi2','xi3','xi4'};
vart = [r_t,xi_t];

for i = 1:6
    figure(2+i)
    boxplot(eval([vars{i},'(:,inds)']),inds);
    line([0,10000],[vart(i),vart(i)],'LineStyle','-.','Color','c');
    xtickangle(45);
    DrawPDF(['Simu_case_1_',vars{i}]);
end

clear
close all
load('MALG_GMM_2C.mat');
figure(1)
rs = [r1(1,:);r2(1,:)]';
r_t = [1/3,2/3];
es =  (rs-r_t);
[S1,AX1,BigAx1,H1,HAx1] = plotmatrix(es);
DrawPNG('MALG_case_1_PlotMatrix_r');

figure(2)
Ss = [s1(1,:);s2(1,:)]';
S_t = [1,4];
Ss = (Ss-S_t);
[S2,AX2,BigAx2,H2,HAx2] = plotmatrix(Ss);
DrawPNG('MALG_case_1_PlotMatrix_s');

figure(3)
xis = [xi1(1,:);xi2(1,:);xi3(1,:);xi4(1,:)]';
xi_t = [8/3,4/3,1/3,2/3];
xis = (xis-xi_t);
[S3,AX3,BigAx3,H3,HAx3] = plotmatrix(xis);
DrawPNG('MALG_case_1_PlotMatrix_xi');

inds = (1:20)*500;

vars = {'r1','r2','s1','s2','xi1','xi2','xi3','xi4'};
vart = [r_t,S_t,xi_t];

for i = 1:8
    figure(3+i)
    boxplot(eval([vars{i},'(:,inds)']),inds);
    line([0,10000],[vart(i),vart(i)],'LineStyle','-.','Color','c');
    xtickangle(45);
    DrawPDF(['MALG_case_1_',vars{i}]);
end


function [ ] = DrawPNG(name)
    a = 1.1;
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 a 0; 0 -1 0 a; 0 0 1 0; 0 0 0 1]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print('-painters','-dpng','-r200',name)
end

function [ ] = DrawPDF(name)
    a = 1.1;
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 a 0; 0 -1 0 a; 0 0 1 0; 0 0 0 1]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print('-painters','-dpdf','-r200',name)
end
