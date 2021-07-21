% Example of algorithms to find geodesics on the space of parameterized curves and 
%  on the space of unparameterized curves 
% authors: Zhengwu Zhang, Bayan Saparbayeva
% emails: zhengwustat@gmail.com, saparbayevabt@gmail.com
% date: Jul. 30, 2021
clear all;
close all;
% add paths
addpath('./Phase-AmplitudeFunction')
addpath('./Phase-Amplitude-modeling')
addpath('./SimulatedData')
% load the data
load SimulationPaths8.mat;
load gammas0.mat
% Get pair of curves p1 and p2 :
p1 = paths0{1};
p2 = paths0{6};
% helps draw 'p1' and 'p2'
T = size(p1, 2); % T=size(p2, 2)
t = linspace(0, 1, T);
% draw the baseline 'beta'
s = linspace(0, 1, 1001);
% tools to draw the geodisics in C
N2 = 18;
s2 = linspace(0, 1, N2); 
% reparameterization
gamma = gammas0{1};
p2 = Group_Action_by_Gamma_p(p2, gamma);
% starting points
x1 = p1(:, 1);
x2 = p2(:, 1);
% TSRVCs
q1 = path_to_q(p1);
q2 = path_to_q(p2);
%
% Extract the set of points for 'theta'
M = 401;
thetaT = linspace(-pi/2, pi/2, M);
for i = 1:M
    SCT2(i) = SquareLength(p1, p2, thetaT(i));
    dSCT2(i) = DthetaSquareLength(p1, p2, thetaT(i));
end
% draw the values of the square length and its derivatives
figure('DefaultAxesFontSize', 26)
plot(thetaT, SCT2, thetaT, dSCT2);
lgd1 = legend({'$d^2(p_0, p, \theta)$',...
 '$\frac{\partial}{\partial\theta}\Big(d^2(p_0, p, \theta)\Big)$'}, ...
    'interpreter','latex');
lgd1.FontSize = 25;
%}
%
% Find the parameter 'theta' for the baseline 'beta' of 
% the geodesic connecting 'p1' and 'p2' on 
% the space of parameterized curves :
[theta1, square_length1, gradient_theta1] =...
    theta_geodesic(p1, p2, 0, 50, 0.0000000000000000001, 60, 30, 30);
%Construct the baseline  on the space of parameterized curves :
normal1 = (x1 + x2)*sin(theta1)/norm(x1 + x2) +...
     cross(x1,x2)*cos(theta1)/norm(cross(x1, x2));
phi1 = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta1)^2-...
             dot(x1, x2)*sin(theta1)^2)));  
beta1 = x1.*cos(s*phi1) + cross(normal1, x1).*sin(s*phi1) + ...
    dot(normal1, x1)*normal1.*(1 - cos(s*phi1));
% Construct the geodesic on the space of parameterized curves 
visparq1 =...
    visparallel_transport_PL(q1, x1, x2, N2, theta1);
visparq2 =...
    visparallel_transport_PL(q2, x2, x1, N2, -theta1);
beta11 = x1.*cos(s2*phi1) + cross(normal1, x1).*sin(s2*phi1) + ...
    dot(normal1, x1)*normal1.*(1 - cos(s2*phi1));
% Plot the baseline between curves 'p1' and 'p2' 
figure(2); clf;
[x, y, z] = sphere(100);
h = surf(0.99*x, 0.99*y, 0.99*z);
axis equal off;
colormap summer;
grid off;
set(h, 'LineStyle', 'none');%, 'FaceAlpha',0.75);
hold on;
pp1 = plot3(p1(1, :), p1(2, :), p1(3, :), ...
     'blue', 'LineWidth', 1.2);
pp2 = plot3(p2(1, :), p2(2, :), p2(3, :), ...
     'cyan', 'LineWidth', 1.2);
bbeta1 = plot3(beta1(1, :), beta1(2, :), beta1(3, :), ...
     'red', 'LineWidth', 1.6);
hold off;
%{
lgd3 = legend([pp1, pp2, bbeta1], ...
   {'$p_0$', '$p$', '$\beta$'}, ...
    'interpreter','latex');
lgd3.FontSize = 21;
 %}
% Draw the geodesics in the space of parameterized curves 
figure(3); clf;
[x, y, z] = sphere(100);
h = surf(0.99*x, 0.99*y, 0.99*z) ;
axis equal off;
colormap summer;
grid off;
set(h, 'LineStyle', 'none');%, 'FaceAlpha',0.75);
hold on;
plot3(beta1(1, :), beta1(2, :), beta1(3, :), ...
    'red', 'LineWidth', 2);
for j=2:(N2-1)
 pgeo(:, :, j) = q_to_path((visparq1(:, :, j) -...
                s2(j)*(visparq1(:, :, j) -...
                visparq2(:, :, N2 - j + 1))), beta11(:, j));
 newpp = pgeo(:, :, j);
 plot3(newpp(1, :), newpp(2, :), newpp(3, :), ...
      'red', 'LineWidth', 1);
end
plot3(p1(1, :), p1(2, :), p1(3, :), ...
 'blue', 'LineWidth', 1.2);
plot3(p2(1, :), p2(2, :), p2(3, :), ...
   'cyan', 'LineWidth', 1.2);
hold off;
%}
% Find the parameter 'theta' for the baseline 'beta'
% of the geodesic connecting 'p1' and 'p2' on 
% the space of unparameterized curves, 
% the reparameterization 'gamman' for 'p1', 
% the optimally reparameterized curve 'p1n',
% the reparameterization 'inv_gamman' for 'p2', 
% the optimally reparameterized curve 'p2n',
% and the square length of the geodesic 'square_length1n' :
[theta1n, p1n, p2n, gamman, inv_gamman, square_length1n, gradient_theta1n] = theta_geodesic_alignment2(p1, p2,...
    0, 50, 5*10^(-8), 30, 15, 15);
% Draw the alignment 'gamman':
figure('DefaultAxesFontSize', 26);clf;
plot(t, gamman);
lgd2 = legend({'$\gamma^{*}$'}, ...
    'interpreter','latex');
hold off;
lgd2.FontSize = 34;
% Construct the baseline 'beta1n' on the space of unparameterized curves :  
normal1n = (x1 + x2)*sin(theta1n)/sqrt(2 + 2*dot(x1, x2)) +...
     cross(x1, x2)*cos(theta1n)/sqrt(1 - dot(x1, x2)^2);
phi1n = 2*asin(sqrt((1 - dot(x1, x2))/...
    (1 + cos(theta1n)^2 - dot(x1, x2)*sin(theta1n)^2)));  
beta1n = x1.*cos(s*phi1n) + cross(normal1n, x1).*sin(s*phi1n) + ...
    dot(normal1n, x1)*normal1n.*(1 - cos(s*phi1n));
% Construct the geodesic on the space of unparameterized curves 
q1n = path_to_q(p1n);
visparq1n =...
    visparallel_transport_PL(q1n, x1, x2, N2, theta1n);
visparq2n =...
    visparallel_transport_PL(q2, x2, x1, N2, -theta1n);
beta11n = x1.*cos(s2*phi1n) + cross(normal1n, x1).*sin(s2*phi1n) + ...
    dot(normal1n, x1)*normal1n.*(1 - cos(s2*phi1n));     
% Draw the geodesics in the space of unparameterized curves
figure(5); clf;
[x, y, z] = sphere(100);
h = surf(0.99*x, 0.99*y, 0.99*z);
axis equal off;
colormap summer;
grid off;
set(h, 'LineStyle', 'none');%, 'FaceAlpha',0.75);
hold on;
plot3(beta1(1, :), beta1(2, :), beta1(3, :), ...
     'red', 'LineWidth', 1.6, 'LineStyle', '--');
plot3(beta1n(1, :), beta1n(2, :), beta1n(3, :), ...
     'color', [0.4940 0.1840 0.5560], 'LineWidth', 1.6);
 for j=2:(N2-1)
  pgeo(:, :, j) = q_to_path((visparq1n(:, :, j) -...
                s2(j)*(visparq1n(:, :, j) -...
                visparq2n(:, :, N2-j+1))), beta11n(:, j));
  newpp = pgeo(:, :, j);
  plot3(newpp(1, :), newpp(2, :), newpp(3, :), ...
      'color', [0.4940 0.1840 0.5560], 'LineWidth', 1);
 end
plot3(p1(1, :), p1(2, :), p1(3, :), ...
     'blue', 'LineWidth', 1.2);
plot3(p2(1,: ), p2(2, :), p2(3, :), ...
     'cyan', 'LineWidth', 1.2); 
hold off
%
% Draw 'beta1' and 'bata1n'  :
figure(6); clf;
[x, y, z] = sphere(100);
h = surf(0.99*x, 0.99*y, 0.99*z) ;
axis equal off;
colormap summer;
grid off;
set(h, 'LineStyle', 'none');%, 'FaceAlpha',0.75);
hold on;
ppp1 = plot3(p1(1, :), p1(2, :), p1(3, :), ...
     'blue', 'LineWidth', 1.2);
ppp2 = plot3(p2(1, :), p2(2, :), p2(3, :), ...
     'cyan', 'LineWidth', 1.2);
bbbeta1 = plot3(beta1(1, :), beta1(2, :), beta1(3, :), ...
     'red', 'LineWidth', 1.6);
bbbeta1n = plot3(beta1n(1, :), beta1n(2, :), beta1n(3, :), ...
     'color', [0.4940 0.1840 0.5560], 'LineWidth', 1.6);
hold off;
%{
lgd4 = legend([ppp1, ppp2, bbbeta1, bbbeta1n], {'$p_0$', ...
     '$p$', '$\beta$', '$\beta^{*}$'}, ...
     'interpreter','latex');
 lgd4.FontSize = 21;
%}
%}
% Clear memory
clear T t N2 x y z j h s s2 thetaT;
clear tmppath newpp pgeo lgd1 lgd2;
clear alpha1 alpha1n phi1 phi1n normal1 normal1n;
clear beta1 beta1n beta11 beta11n dbeta1n0 dbeta10;
clear p1 p1n p2 q1 q1n q2 x1 x2;
clear e2n e2 e1n e1 i M SCT2 dSCT2;
clear gamma gamman gammas0; 
clear visparq1 visparq1n visparq2 visparq2n;
clear epu1 epw1 theta1n p2n inv_gamman;
clear bbbeta1 bbbeta1n bbeta1 length1 length2Our;
clear pp1 pp2 ppp1 ppp2 theta1;
