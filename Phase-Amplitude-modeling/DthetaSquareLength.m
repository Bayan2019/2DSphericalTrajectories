function [gradient_theta, square_length] = DthetaSquareLength(p1, p2, theta)
%%% the derivative 'gradient_theta' of the square length 'SquareLength'
%%% with respect to 'theta' on the space of 
%%% parameterized curves
%%% Inputs: two curves 'p1' and 'p2', 'theta' which determines 
%%%         the baselines 'beta'
%%% Outputs: the derivative 'gradient_theta' 
%%%          of the square length 'SquareLength'
%%%          with respect to 'theta'
%%%          and the square length 'square_length'
% if we assume that q1 and q2 have the same length 
% or determined at the same finite set of points in [0, 1]
T = size(p1, 2); % T=size(p2, 2)
% we have starting points on S2
x1 = p1(:, 1);
x2 = p2(:, 1);
% and SRVC at tangent spaces
q1 = path_to_q(p1);
q2 = path_to_q(p2);
% 'alpha' and 'phi'
alpha = acos(sin(theta)*sqrt((1 + dot(x1, x2))/2));
phi = 2*asin(sqrt((1 - dot(x1, x2))/(1 + cos(theta)^2 -...
             dot(x1, x2)*sin(theta)^2)));
% parallel transport of 'q1' from 'x1' to 'x2'
parq1 = parallel_transport_PL(q1, x1, x2, theta);
% computing square of length         
len2 = phi^2*sin(alpha)^2;  
SC = trapz(linspace(0, 1, T), sum((parq1 - q2).*(parq1 - q2)));
square_length = len2 + SC;
% computing the gradient
Xix = dot(x1, x2); 
dalpha_theta = -(cos(theta)*(Xix/2 + 1/2)^(1/2))/...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2); 
dphi_theta = -(sin(2*theta)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2))/...
    ((1 - sin(theta)^2)^(1/2)*(sin(theta)^2 + Xix*sin(theta)^2 - 2));
gradient1 = 2*phi*sin(alpha)*...
    (dphi_theta*sin(alpha) + phi*dalpha_theta*cos(alpha));
parq1_theta = d_parallel_transport_PL_theta(q1, x1, x2, theta);
gradient_theta = gradient1 +...
    2*trapz(linspace(0, 1, T), sum((parq1_theta).*(parq1 - q2)));
end