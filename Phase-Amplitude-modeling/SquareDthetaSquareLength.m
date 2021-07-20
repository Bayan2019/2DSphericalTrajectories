function gradient_theta2 = SquareDthetaSquareLength(p0, u, w, vartheta)
%%% 
%%% Inputs: the curve 'p0', the direction 'u' 
%%%         of the baseline 'beta' which determines 
%%%         the differece 'w' of TSRVCs,
%%%         the value 'vartheta' that determines 
%%%         the baseline 'beta' and also the second curve 'p'
%%%         with the help of the direction 'u' and TSRVC difference 
%%% Outputs: 'gradient_theta2' is the square  of 'DthetaSquareLength'
 x0 = p0(:, 1);
 q0 = path_to_q(p0);
 r = norm(u);
 normal = x0*cos(vartheta) + ...
     (cross(x0, u)/norm(cross(x0, u)))*sin(vartheta);
 phi = r/sin(vartheta);
 x = x0*cos(phi) + cross(normal, x0)*sin(phi) +...
     normal*dot(normal, x0)*(1 - cos(phi));
 theta = asin(sqrt(2/(1 + dot(x0, x)))*cos(vartheta));
 q = parallel_transport_PL(q0+w, x0, x, theta);
 p = q_to_path(q, x);
 d_theta_square_length = DthetaSquareLength(p0, p, theta);
 gradient_theta2 = d_theta_square_length^2;
end

