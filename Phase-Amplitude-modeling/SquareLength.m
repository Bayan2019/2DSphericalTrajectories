function square_length = SquareLength(p1, p2, theta)
%%% the square length 'SquareLength' 
%%% on the space of parameterized curves
%%% Inputs: two curves 'p1' and 'p2', 'theta' which determines 
%%%         the baselines 'beta'
%%% Outputs: the square of length 'square_length' on
%%%          the space of parameterized curves
% if we assume that q1 and q2 have the same length 
% or determined at the same finite set of points in [0, 1]
 T = size(p1, 2); % T=size(p2, 2)
% we have starting points on S2
x1 = p1(:, 1);
x2 = p2(:, 1);
% and SRVC at tangent spaces
q1 = path_to_q(p1);
q2 = path_to_q(p2);
alpha = acos(sin(theta)*sqrt((1+dot(x1, x2))/2));
phi = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta)^2-...
             dot(x1, x2)*sin(theta)^2)));
len2 = phi^2*sin(alpha)^2;
%SC = 0;
parq1 = parallel_transport_PL(q1, x1, x2, theta);
SC = trapz(linspace(0,1,T), sum((parq1-q2).*(parq1-q2)));
square_length = SC + len2;
end

