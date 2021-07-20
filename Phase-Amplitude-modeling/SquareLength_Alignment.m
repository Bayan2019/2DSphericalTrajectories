function [square_length, p1n, p2n, gamma, inv_gamma] = SquareLength_Alignment(p1, p2, theta)
%%% the square length 'SquareLength_Alignment' 
%%% on the space of unparameterized curves
%%% Inputs: two curves 'p1' and 'p2', 'theta' which determines 
%%%         the baselines 'beta'
%%% Outputs: the square of length 'square_length' on the space of
%%%          unparameterized curves
%%%          the reparameterization 'gamman' for 'p1', 
%%%          the optimally reparameterized curve 'p1n',
%%%          the reparameterization 'inv_gamman' for 'p2', 
%%%          the optimally reparameterized curve 'p2n',
% if we assume that q1 and q2 have the same length 
% or determined at the same finite set of points in [0, 1]
 %T = size(p1, 2); % T=size(p2, 2) 
% we have starting points on S2
x1 = p1(:, 1);
x2 = p2(:, 1);
% and SRVC at tangent spaces
q1 = path_to_q(p1);
q2 = path_to_q(p2);
% parallel transport
par_q1 = parallel_transport_PL(q1, x1, x2, theta);
% alignment of the parallel transport parq1
T = size(p1, 2); % T=size(p2, 2)
t = linspace(0,1,T);
[G,T1] = DynamicProgrammingQ2(q2, t, par_q1, t, t, t, 0);
gamma = interp1(T1, G, t);
p1n = Group_Action_by_Gamma_p(p1, gamma);
inv_gamma = invertGamma(gamma);
inv_gamma = (inv_gamma-inv_gamma(1))/(inv_gamma(end)-inv_gamma(1));
p2n = Group_Action_by_Gamma_p(p2, inv_gamma);
% computing the length
square_length1n = SquareLength(p1n, p2, theta);
square_length2n = SquareLength(p1, p2n, theta);
square_length = min(square_length1n, square_length2n);
end