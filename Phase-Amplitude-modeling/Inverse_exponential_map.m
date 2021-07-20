function [u, w, theta, square_length, gradient_theta] = Inverse_exponential_map(p0, p, ...
    theta0, step, epsilon, NN1, NN2, NN3)
%%% the Inverse Exponential map on the space of parameterized curves
%%% the function is finding 'theta' that determines the baseline which generates
%%% the geodesic between p1 and p2 and the space of parameterized curves
%%% the inputs: two curves 'p1' and 'p2', the initial parameter '\theta',
%%%             the degree of accurace 'epsilon', the step size 'step',
%%%             the maximum number of iterations 'NN1',
%%%             the maximum number of iterations when we don't perform
%%%             steps 'NN2' and 'NN3'
%%% the outputs: the 'theta' that determines baseline 'beta',
%%%              the optimal square length 'square_length',
%%%              the gradient of square length 'gradient_theta' 
%%%              with respect to 'theta',
%%%              the '(u, w)', the direction 'u' of the baseline 'beta'
%%%              the differece 'w' of TSRVCs
% find optimal 'theta' which determines 'beta'
 [theta, square_length, gradient_theta] = theta_geodesic(p0, p, ...
     theta0, step, epsilon, NN1, NN2, NN3);
% starting points:
 x0 = p0(:, 1);
 x = p(:, 1);
 phi = 2*asin(sqrt((1 - dot(x0, x))/(1 + cos(theta)^2 -...
             dot(x0, x)*sin(theta)^2)));
 normal = (x0 + x)*sin(theta)/norm(x0 + x) +...
     cross(x0, x)*cos(theta)/norm(cross(x0, x));
% direction of the p-optimal curve at the beginning
 u = -x0*phi*sin(0) + cross(normal, x0)*phi*cos(0) +...
       normal*dot(normal, x0)*phi*sin(0);
 q0 = path_to_q(p0);
 q = path_to_q(p);  
 parq = parallel_transport_PL(q, x, x0, -theta);
% the difference 
 w = parq - q0;
end

