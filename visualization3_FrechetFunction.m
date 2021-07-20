% Example of implementing algorithm to find the sample Frechet mean of
% the set of curves (simulated)
% on the space of unparameterized curves 
% authors: Zhengwu Zhang, Bayan Saparbayeva
% emails: zhengwustat@gmail.com, saparbayevabt@gmail.com
% date: Jul. 30, 2021
clear all;
close all;
% add paths
addpath('./amplitude_separation')
addpath('./Bayan/denseFD_for_unparameterized_curves')
addpath('./Bayan/denseFD_for_unparameterized_curves/SimulatedCurves/5')
% load curves 
load SimulatedCurves2.mat;
%
[mup,muq,mupath, KF_30_60] = KarcherMean(paths0, 'slow', 30, 60);
[mup,muq,mupath, KF_30_120] = KarcherMean(paths0, 'slow', 30, 120);
[mup,muq,mupath, KF_30_240] = KarcherMean(paths0, 'slow', 30, 240);
[mup,muq,mupath, KF_60_60] = KarcherMean(paths0, 'slow', 60, 60);
[mup,muq,mupath, KF_60_120] = KarcherMean(paths0, 'slow', 60, 120);
[mup,muq,mupath, KF_60_240] = KarcherMean(paths0, 'slow', 60, 240);
[mup,muq,mupath, KF_120_60] = KarcherMean(paths0, 'slow', 120, 60);
[mup,muq,mupath, KF_120_120] = KarcherMean(paths0, 'slow', 120, 120);
[mup,muq,mupath, KF_120_240] = KarcherMean(paths0, 'slow', 120, 240);
%}
% Find the sample Frechet mean 'p' with the starting point 'x0'
% vector 'theta' of parameters determining the baselines 'beta's
% connecting 'p' and 'pj' from 'paths0'
% set of curves 'paths' that are aligned curves from 'paths0' to 'p' 
% on the space of unparameterized curves, 
% the reparameterization 'gammas' of 'paths0' to 'paths',
% the reparameterization 'inv_gammas' of 'paths' to 'paths0',
% and the optimal sample Frechet function 'SFF':
[p, x0, theta, paths, gammas, inv_gammas, SFF] = SFrechet_Alignment2(paths0, ...
    30000.6, 50, 0.0001, 0.0000001, 0.00001, 60, 30, 30);
n = 6;
disp(num2str(SFF, n));
%
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
%start drawing
% T = size(p, 2); % T=size(p2, 2)
N = length(paths0);
% Draw curves:
figure(1)
[x, y, z] = sphere(100);
h = surf(0.99*x, 0.99*y, 0.99*z) ;
axis equal off;
colormap summer;
grid off;
set(h, 'LineStyle', 'none');
hold on;
for j = 1:N
 tmppath = paths0{j};
 plot3(tmppath(1, :), tmppath(2, :), tmppath(3, :), 'blue', 'linewidth', 1.0);
end
%title('Original Samples');
hold off;
%}
% Clear memory
clear x y z j h N;
clear tmppath n;
clear x0 theta paths inv_gammas0 gammas0;
clear mup muq mupath;

