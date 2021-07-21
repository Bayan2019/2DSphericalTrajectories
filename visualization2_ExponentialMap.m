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
% reparameterization
gamma = gammas0{1};                   
p2 = Group_Action_by_Gamma_p(p2, gamma);
% implement inverse exponential map
[u, w, theta1, square_length, gradient_theta] = Inverse_exponential_map(p1, p2, 0, 50, 0.000000000000001, 60, 30, 30);
% implemen the exponential map
[pp2, theta2, vartheta2, gradient_theta2] = exponential_map(p1, u, w, 50, 0.0000000000001, 60, 30, 30);
% compare the outcomes
L = SquareLength(p2, pp2, 0);
% clear the memory
clear gamma gammas0 gradient_theta gradient_theta2 vartheta2 p1 p2;
clear paths0 pp2 t T theta1 theta2 u w square_length;
