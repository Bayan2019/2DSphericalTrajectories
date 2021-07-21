% Example of implementing the exponential map and 
%     the inverse exponential on the space of parameterized curves 
% authors: Zhengwu Zhang, Bayan Saparbayeva
% emails: zhengwustat@gmail.com, saparbayevabt@gmail.com
% date: Jul. 30, 2021
clear all;
close all;
% add paths
addpath('./Phase-AmplitudeFunction')
addpath('./Phase-Amplitude-modeling')
addpath('./SimulatedData')
% load data
load SimulationPaths8.mat;
load gammas0.mat
% Get pair of curves p1 and p2 :
p1 = paths0{11};
p2 = paths0{51};
% Starting points
x1 = p1(:, 1);
x2 = p2(:, 1);
% reparameterization
gamma = gammas0{1};                   
p2 = Group_Action_by_Gamma_p(p2, gamma);
% implement inverse exponential map
[u1, w1, theta1, square_length, gradient_theta] = Inverse_exponential_map(p1, p2, 0, 50, 0.000000000000001, 60, 30, 30);
% implemen the exponential map
[pp2, theta2, vartheta2, gradient_theta2] = exponential_map(p1, u1, w1, 50, 0.0000000000001, 60, 30, 30);
% compare the outcomes
L = SquareLength(p2, pp2, 0);
% basis in the tangent plane at x1        
 e_x1 = cross(x1, x2)/norm(cross(x1, x2));
 e_x2 = cross(cross(x1, x2), x1)/norm(cross(cross(x1, x2), x1));
% Define grids on the tangent space at x1
 [u, v] = meshgrid(-2.5:1/80:2.5);
% the tangent space at x1
 tangentplane1 = @(u, v) u*e_x1(1) + v*e_x2(1) + x1(1);
 tangentplane2 = @(u, v) u*e_x1(2) + v*e_x2(2) + x1(2);
 tangentplane3 = @(u, v) u*e_x1(3) + v*e_x2(3) + x1(3); 
% projection of 'u1' and 'w1' to the tangent space at x1 (which are exactly 'u1' and 'w1')
 u11 = e_x1'*(t.*u);
 u12 = e_x2'*(t.*u);
 ww11 = e_x1'*w;
 ww12 = e_x2'*w;
% plot the image '(u, w) = exp_{p_1}^{-1}p2' of the inverse exponential map at the tangent space at 'p1'
 figure('DefaultAxesFontSize', 26); clf;
 plot(u11, u12, 'color', [0.3010 0.7450 0.9330], 'LineStyle', '--', ...
     'LineWidth', 1.8)
 hold on;
 uc = plot(u11(T), u12(T), 'w', 'Marker', '*', ...
     'MarkerEdgeColor', [0.3010 0.7450 0.9330], 'LineWidth', 1.9, ...
     'MarkerSize', 9.4);
 www = plot(w11, w12, 'color', [0.3010 0.7450 0.9330], 'LineWidth', 1.8);
 lgd3 = legend([uc, www],...
   {'$u$', '$w$'}, 'interpreter','latex');
 lgd3.FontSize = 35;
 hold off;
% connect 'u1' and 'w1' to the tangent space at x1
 ww = w1 + x1;
 uu = t.*u1 + x1;
% plot '(u1, w1')' at the tangent spact at x1
 figure(5); clf;
 [x, y, z] = sphere(100);
 h = surf(0.999*x, 0.999*y, 0.999*z) ;
 axis equal off;
 colormap summer;
 grid off;
 set(h, 'LineStyle','none');
 hold on;
 fsurf(tangentplane1, tangentplane2, tangentplane3,...
     [-0.6 1.6 -0.6, 0.5], 'FaceAlpha', 0.46, ...
     'EdgeColor','none', 'FaceColor', 'blue');
 plot3(p1(1, :), p1(2, :), p1(3, :), 'blue', 'LineWidth', 1.2);
 plot3(x1(1), x1(2), x1(3),'b*', 'LineWidth', 1.6);
 plot3(uu(1, T), uu(2, T), uu(3, T),'c*', 'LineWidth', 9.6);
 plot3(uu(1, :), uu(2, :), uu(3, :), 'cyan', 'LineWidth', 5.5, 'LineStyle', '--');
 plot3(ww(1,:), ww(2,:), ww(3,:), 'cyan', 'LineWidth', 1.5);
 hold off;
% clear the memory
clear e_x1 e_x2 x1 x2
clear gamma gammas0 gradient_theta gradient_theta2 vartheta2 p1 p2;
clear paths0 pp2 t T theta1 theta2 u w square_length;
