function SFF = SampleFrechetFunction(x, theta, paths)
%%% the sample Frechet function on the space of parameterized curves
%%% Inputs: 'x' is the starting point of the sample Frechet mean,
%%%         'paths' is the set of curves, 
%%%         'theta' is vector of parameters 'thetaj' that determine 
%%%            the geodesic baselines connecting 'x' with 'xj',         
%%% Outputs: the sample Frechet fucntion 'SFF'
% sample size
N = length(theta);
% if we assume that q1 and q2 have the same length 
% or determined at the same finite set of points in [0, 1]
%[d,T] = size(patharray{1}); % T=size(p2, 2)
% we have starting points on S2
q=0;
for j=1:N
pj = paths{j};
Xj = pj(:, 1);
qj = path_to_q(pj);
parqj = parallel_transport_PL(qj, Xj, x, theta(j));
q = q + parqj;
end
q = q/N;
p = q_to_path(q, x);
SFF = 0;
for j=1:N
pj = paths{j};
SFF = SFF + SquareLength(pj, p, theta(j));
end
SFF = SFF/N;
end