function [gradient_x, SFF] = DxSampleFrechetFunction(x, theta, paths)
%%% the derivative 'gradient_x' of the sample Frechet function 
%%% 'SampleFrechetFunction' with respect to 'x' on the space of 
%%% parameterized curves
%%% Inputs: 'x' is the starting point of the sample Frechet mean,
%%%         'paths' is the set of curves, 
%%%         'theta' is vector of parameters 'thetaj' that determine 
%%%            the geodesic baselines connecting 'x' with 'xj',         
%%% Outputs: the derivative 'gradient_x' of 'SampleFrechetFunction'
%%%          with respect to 'x'
%%%          and the sample Frechet fucntion 'SFF'
% sample size
N = length(theta);
% if we assume that qis have the same length 
% or determined at the same finite set of points in [0, 1]
%[d,T] = size(patharray{1}); % T=size(p2, 2)
% we have starting points on S2
q = 0;
for j = 1:N
pj = paths{j};
Xj = pj(:, 1);
qj = path_to_q(pj);
parqj = parallel_transport_PL(qj, Xj, x, theta(j));
q = q + parqj;
end
q = q/N;
p = q_to_path(q, x);
gradient_x = [0; 0; 0];
SFF = 0;
for j=1:N
pj = paths{j};
[gradient_xj, square_lengthj] = DxSquareLength(pj, p, theta(j));
gradient_x = gradient_x + gradient_xj;
SFF = SFF + square_lengthj;
end
gradient_x = gradient_x/N;
SFF = SFF/2;
end