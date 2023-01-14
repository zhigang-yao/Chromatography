% Compare 3 samplers
clear
clc
close all

load data.mat
[xi_MH,~] = MMH(y_obs, Tn, 10, 10000);
[xi_MGDG,~,~,~]=MGDG(y_obs, Tn, 10, 10000,'Loss-ratio.mat');
[xi_MALG,~,~,~]=MGDG(y_obs, Tn, 10, 10000,'Loss-ratio.mat');

figure(1)
plot(sort(xi_MGDG,2))
ylim([0,3.5])
DrawPDF('Simu_case_1_MGDG')
figure(2)
plot(sort(xi_MALG,2))
ylim([0,3.5])
DrawPDF('Simu_case_1_MALG')
figure(3)
plot(sort(xi_MH,2))
ylim([0,3.5])
DrawPDF('Simu_case_1_MH')

function [ ] = DrawPDF(name)
    a = 1.1;
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 a 0; 0 -1 0 a; 0 0 1 0; 0 0 0 1]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print('-painters','-dpdf','-r200',name)
end
