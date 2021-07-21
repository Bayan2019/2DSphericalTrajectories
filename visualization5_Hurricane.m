% Example of calculating sample Frechet mean on the space of parameterized curves and 
%  on the space of unparameterized curves for Hurricane data
% authors: Zhengwu Zhang, Bayan Saparbayeva
% emails: zhengwustat@gmail.com, saparbayevabt@gmail.com
clear all;
close all;
% add paths
addpath('./Phase-AmplitudeFunction')
addpath('./Phase-Amplitude-modeling')
addpath('./RealData')
% load data
load ProcessHurricanePaths.mat;
% variables for curves
T = 100;
% time grid
t = linspace(0, 1, T);
N = 10; % the number of samples
for i=1:N
 Resampled_track = ReSampleSphereTraj(PPProcess_Hurri_Path{2*i - 1}, T);
 Smoothed_track = SmoothPath(Resampled_track, 7, 1);
 paths0{i} = Smoothed_track;
end
% calculating the sample Frechet function on the space of parameterized
% curve
[p1, x1, theta1] = SFrechet(paths0, 30000.6, 50, 10^(-7), 10^(-6), 60, 30, 30);
% calculating the sample Frechet function on the space of unparameterized
% curves
[p2, x2, theta2, paths2, gammas, inv_gammas, FF] = SFrechet_Alignment2(paths0, 30000.6, ...
    50, 0.0001, 10^(-7), 10^(-5), 60, 30, 30);  
% variables for variances
for j =1:N
 pj = paths0{j};
 [uj, wj, thetaj, SCj, gradientj] = Inverse_exponential_map(p1, pj, -theta1(j), 50, 10^(-9), 60, 25, 25);
 u1{j} = uj;
 w1{j} = wj;
end
% variables for aligned variances
for j =1:N
 pj = paths2{j};
 [uj, wj, thetaj, SCj, gradientj] = Inverse_exponential_map_Alignment(p2, pj, -theta2(j), 50, 10^(-8), 30, 15, 15);
 u2{j} = uj;
 w2{j} = wj;
end 
% draw the sample Frechet mean on the space of parameterized curves
figure(1);clf;
axis equal off;
grid off;
hold on;
globe([], 'earth_1600.png');
hold on;
for j=1:N
 pj = paths0{j};
 ppj = plot3(pj(1, :), pj(2, :), pj(3, :), 'yellow', 'linewidth', 0.9);
end
pmu = plot3(p1(1, :), p1(2, :), p1(3, :), 'color', [1 0.6 0], 'LineWidth', 3);
lgd1 = legend([pmu, ppj], {'$p_{\mu}$', '$p_j$'}, ...
   'interpreter','latex');
lgd1.FontSize = 21;
%title('Frechet mean in $\mathbb{C}$', 'interpreter','latex');
hold off;
% drawing the sample Frechet mean on the space of unparameterized curves 
figure(2)
axis equal off;
grid off;
hold on;
globe([], 'earth_1600.png');
%arg='100, ''yellow'', ''fill'', ''markeredgecolor'', ''black'' ';
hold on;
for j=1:N
 pj = paths0{j};
 ppj = plot3(pj(1, :), pj(2, :), pj(3, :), 'yellow', 'linewidth', 0.9);
end
pmu = plot3(p2(1, :), p2(2, :), p2(3, :), 'color', [0.9 0.3 1], 'LineWidth', 3);
lgd2 = legend([pmu, ppj], {'$\tilde{p}_{\mu}$', '$p_j$'}, ...
   'interpreter','latex');
lgd2.FontSize = 21;
%title('The Frechet mean in $\mathbb{C}/\Gamma$', 'interpreter','latex');
hold off;
% drawing reparameterizations
for i=1:N
    Gammas0(:, i) = gammas{i};
end
figure('DefaultAxesFontSize', 26); clf;
plot(t, Gammas0);
title('Phase variability', 'interpreter', 'latex');
hold off;
% drawing variances
figure(103); clf;
axis equal off;
grid off;
hold on;
globe([], 'earth_1600.png');
%arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
plot3(p1(1, :), p1(2, :), p1(3, :), 'color', [1 0.6 0], 'LineWidth', 3);
%title('Karcher mean and variance along the mean')
for k=1:T   
 for j=1:N
  %uj = u1{j};   
  wj = w1{j};
  Dat1(j, :) = wj(:, k);
 end
 K1 = cov(Dat1);
 [U1, S1, V1] = svd(K1);
 Total_var1(k) = trace(K1);
 if mod(k, 15)==0
  lam1 = sqrt(S1(1, 1));
  lam2 = sqrt(S1(2, 2));
  the = linspace(0, 2*pi, 100);
  [xthe] = 0.6*lam1*sin(the);
  [ythe] = 0.6*lam2*cos(the);
  yyy = U1*[xthe; ythe; zeros(1, 100)];% + repmat(p1(:, k), 1, 100);
  for i = 1:100
  yyy(:, i) = ForwardParallelTranslation(p1(:, 1:k), yyy(:, i)) + p1(:, k); 
  end  
  plot3(yyy(1, :), yyy(2, :), yyy(3, :), 'y', 'LineWidth', 3);
 end
end
hold off;
% drawing aligned variances
figure(104); clf;
axis equal off;
grid off;
hold on;
globe([], 'earth_1600.png');
hold on;
plot3(p2(1, :), p2(2, :), p2(3, :), 'color', [0.9 0.3 1], 'LineWidth', 3);
%title('Karcher mean and variance along the mean')
for k=1:T   
 for j=1:N
  %uj = u2{j};   
  wj = w2{j};
  Dat2(j, :) = wj(:, k);
 end
 K2 = cov(Dat2);
 [U2, S2, V2] = svd(K2);
 Total_var_align2(k) = trace(K2);
 if mod(k, 15)==0
  lam1 = sqrt(S2(1, 1));
  lam2 = sqrt(S2(2, 2));
  the = linspace(0, 2*pi, 100);
  [xthe] = 0.6*lam1*sin(the);
  [ythe] = 0.6*lam2*cos(the);
  yyy = U2*[xthe; ythe; zeros(1, 100)];% + repmat(p2(:, k), 1, 100);
  for i = 1:100
  yyy(:, i) = ForwardParallelTranslation(p2(:, 1:k), yyy(:, i)) + p2(:, k);
  end     
  plot3(yyy(1, :), yyy(2, :), yyy(3, :), 'y', 'LineWidth', 3);
 end
end
hold off;
% draw to compare variances
figure('DefaultAxesFontSize', 26); clf;
pv1 = plot(t, Total_var1, 'color', [1 0.6 0], 'linewidth', 2.5, 'LineStyle', '--');
hold on;
pv2 = plot(t, Total_var_align2, 'color', [0.9 0.3 1], 'linewidth',2.5, 'LineStyle', '--');
lgd3 = legend([pv1, pv2], {'$\rho_{\mu}$', '$\tilde{\rho}_{\mu}$'}, 'interpreter', 'latex');
lgd3.FontSize = 21;
hold off;
%{
% comparison the sample Frechet function
[mup,muq,mupath, KF_30_60] = KarcherMean(paths0, 'slow', 30, 60);
[mup,muq,mupath, KF_30_120] = KarcherMean(paths0, 'slow', 30, 120);
[mup,muq,mupath, KF_30_240] = KarcherMean(paths0, 'slow', 30, 240);
[mup,muq,mupath, KF_60_60] = KarcherMean(paths0, 'slow', 60, 60);
[mup,muq,mupath, KF_60_120] = KarcherMean(paths0, 'slow', 60, 120);
[mup,muq,mupath, KF_60_240] = KarcherMean(paths0, 'slow', 60, 240);
[mup,muq,mupath, KF_120_60] = KarcherMean(paths0, 'slow', 120, 60);
[mup,muq,mupath, KF_120_120] = KarcherMean(paths0, 'slow', 120, 120);
[mup,muq,mupath, KF_120_240] = KarcherMean(paths0, 'slow', 120, 240);
n = 6;
disp(num2str(FF, n));
disp(num2str(KF_30_60, n));
disp(num2str(KF_30_120, n));
disp(num2str(KF_30_240, n));
disp(num2str(KF_60_60, n));
disp(num2str(KF_60_120, n));
disp(num2str(KF_60_240, n));
disp(num2str(KF_120_60, n));
disp(num2str(KF_120_120, n));
disp(num2str(KF_120_240, n));
%}

clear T N t i j Gammas0;
clear uj wj pj thetaj SCj gradientj;
clear arg color ppj ppn pv1 pv2 lgd1 lgd2 lgd3;
clear k Dat1 K1 U1 S1 V1 Total_var1;
clear the xthe ythe the lam1 lam2 yyy;
clear Dat2 K2 U2 S2 V2 Total_var_align2 paths2 paths0;
clear Resampled_track Smoothed_track;
clear theta1 theta2 gammas inv_gammas;
clear mupath muq mup n;
clear KF_30_60 KF_30_120 KF_30_240;
clear KF_60_60 KF_60_120 KF_60_240;
clear KF_120_240 KF_120_120 KF_120_60;
