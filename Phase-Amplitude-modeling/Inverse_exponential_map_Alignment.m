function [u, w, theta, square_length, gradient] = Inverse_exponential_map_Alignment(p0, p, theta0, step, epsilon, NN1, NN2, NN3)
 [theta, p0n, pn, gamma, inv_gamma, square_length, gradient] = theta_geodesic_alignment2(p0, p, ...
   theta0, step, epsilon, NN1, NN2, NN3);
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
 qn = path_to_q(pn);  
 parq = parallel_transport_PL(qn, x, x0, -theta);
 w = parq - q0;
end