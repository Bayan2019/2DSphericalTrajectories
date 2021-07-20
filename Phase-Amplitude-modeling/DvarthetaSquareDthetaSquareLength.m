function [gradient_vartheta, gradient_theta2] = DvarthetaSquareDthetaSquareLength(p0, u, w, vartheta)
%%% 
%%% Inputs: the curve 'p0', the direction 'u' 
%%%         of the baseline 'beta' which determines 
%%%         the differece 'w' of TSRVCs,
%%%         the value 'vartheta' that determines 
%%%         the baseline 'beta' and also the second curve 'p'
%%%         with the help of the direction 'u' and TSRVC difference 
%%% Outputs: the derivative 'gradient_vartheta'
%%%          of square of 'DthetaSquareLength' 
%%%          with respect to 'vartheta' 
%%%          which determines baseline 'beta' and also second curve 'p';
%%%          'gradient_theta2' is the square  of 'DthetaSquareLength' 
 T = size(p0, 2);
 % starting point of 'p0'
 x0 = p0(:, 1);
 % TSRVC of 'p0'
 q0 = path_to_q(p0);
 % the length of 'beta'
 r = norm(u);
 % 'alpha' and 'phi'
 phi = r/sin(vartheta);
 alpha = vartheta;
 % 'normal'
 normal = x0*cos(vartheta) +...
     (cross(x0, u)/norm(cross(x0, u)))*sin(vartheta);
 % starting point of 'p'
 x = x0*cos(phi) + cross(normal, x0)*sin(phi) +...
     normal*dot(normal, x0)*(1 - cos(phi));
 % the angle 'theta' that determines 'beta'
 theta = asin(sqrt(2/(1 + dot(x0, x)))*cos(vartheta));
 % derivatives of 'alpha' and 'phi' with respect 'vartheta'
 dphi_vartheta = (r*cos(vartheta))/(cos(vartheta)^2 - 1);
 dalpha_vartheta = 1;
  % derivative of 'normal' with respect to 'vartheta'
 dnormal_vartheta = -x0*sin(vartheta) +...
     (cross(x0, u)/norm(cross(x0, u)))*cos(vartheta);
 % first derivative of 'alpha' annd 'phi'
 % with respect to 'theta'
 Xix = dot(x0, x);  
 dalpha_theta = -(cos(theta)*(Xix/2 + 1/2)^(1/2))/...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2); 
 dphi_theta = -(sin(2*theta)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2))/...
    ((1 - sin(theta)^2)^(1/2)*(sin(theta)^2 + Xix*sin(theta)^2 - 2));
 x_vartheta = -x0*sin(phi)*dphi_vartheta +...% the second term
     cross(dnormal_vartheta, x0)*sin(phi) +...
     cross(normal, x0)*cos(phi)*dphi_vartheta +... % the third term
     dnormal_vartheta*dot(normal, x0)*(1 - cos(phi)) +...
     normal*dot(dnormal_vartheta, x0)*(1 - cos(phi)) +...
     normal*dot(normal, x0)*sin(phi)*dphi_vartheta;
 % the angle 'theta' that determines 'beta'
 % theta = asin(sqrt(2/(1 + dot(x0, x)))*cos(btheta));
 % the partial derivatives of 'theta' with respect to 'vartheta'
 % and 'x0x'
 x0x = Xix;
 theta_Vartheta = -(2^(1/2)*sin(vartheta)*(1/(x0x + 1))^(1/2))/...
     (1 - (2*cos(vartheta)^2)/(x0x + 1))^(1/2);
 theta_x0x = -(2^(1/2)*cos(vartheta))/(2*(1 - (2*cos(vartheta)^2)/...
     (x0x + 1))^(1/2)*(1/(x0x + 1))^(1/2)*(x0x + 1)^2);
 % final derivative of 'theta' with respect to 'vartheta'
 theta_vartheta = theta_Vartheta + theta_x0x*dot(x0, x_vartheta);
 % second derivatives of 'alpha' and 'phi'
 % with respect to 'theta'
 d2_alpha_theta = -(2*sin(theta)*(Xix/2 - 1/2)*(Xix + 1)^(1/2))/...
    (- sin(theta)^2 - Xix*sin(theta)^2 + 2)^(3/2); 
 d2_phi_theta = (2*(1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*...
    (Xix + 1)^(1/2)*(sin(theta)^2 + Xix*sin(theta)^2 + 2))/...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2)^2; 
 % second derivative of 'alpha' and 'phi' 
 % with respect to 'x' along the direction 'x_vartheta' and 'theta'
 d2_alpha_xutheta = ((cos(theta)*(sin(theta)^2*(Xix/2 + 1/2) - 1) -...
    cos(theta)*sin(theta)^2*(Xix/2 + 1/2))/...
    (4*(Xix/2 + 1/2)^(1/2)*...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(3/2)))*dot(x0, x_vartheta); 
 d2_phi_xutheta = ((cos(theta)*sin(theta)*(2*Xix + 2)*...
    (sin(theta)^2 -...
    2*Xix + Xix*sin(theta)^2))/...
    ((1 - sin(theta)^2)^(1/2)*...
    (1 - Xix)^(1/2)*(Xix + 1)^(3/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2)^2))*dot(x0, x_vartheta);
 % TSRVC of 'p' 
 q = parallel_transport_PL(q0 + w, x0, x, theta);
 % p = q_to_path(q, x);
 % d_theta_square_length = DthetaSquareLength(p0, p, theta);
 % gradient_theta2 = d_theta_square_length^2 
 % derivative of square length from 'x0' to 'x' 
 % with respect to 'theta'
 %
 % parallel transport of 'q0' from 'x0' to 'x'
 parq0 = parallel_transport_PL(q0, x0, x, theta);
 % derivative of 'parq0' with respect to 'theta'
 parq0_theta = d_parallel_transport_PL_theta(q0, x0, x, theta);
 % TSRVC of 'p' 
 % q = parallel_transport_PL(q0 + w, x0, x, theta);
 q_vartheta = d_parallel_transport_PL_xu(q0 + w, x0, x, theta, x_vartheta) +...
     d_parallel_transport_PL_theta(q0 + w, x0, x, theta)*theta_vartheta;
 % computing 'square_gradient_theta'
 d_theta_square_length1 = 2*phi*sin(alpha)*...
    (dphi_theta*sin(alpha) + phi*dalpha_theta*cos(alpha));
 d_theta_square_length2 = 2*...
     trapz(linspace(0, 1, T), sum((parq0_theta).*(parq0 - q)));
 d_theta_square_length = d_theta_square_length1 +...
   d_theta_square_length2;
 gradient_theta2 = d_theta_square_length^2;
 % first derivative of 'dalpha_theta' and 'dphi_theta'
 % with respect to 'vartheta'
 d_dalpha_theta_vartheta = d2_alpha_theta*theta_vartheta +...
    d2_alpha_xutheta;
 d_phi_theta_vartheta = d2_phi_theta*theta_vartheta +...
    d2_phi_xutheta;
 % derivatives of 'd_theta_square_length1' with respect to 
 % 'alpha', 'phi', 'dalpha_theta', and 'dphi_theta'
 d_d_theta_square_length1_alpha = 2*phi*(dphi_theta*sin(2*alpha) +...
    dalpha_theta*phi*cos(2*alpha));
 d_d_theta_square_length1_phi = 2*sin(alpha)*(dphi_theta*sin(alpha) +...
    2*dalpha_theta*phi*cos(alpha));
 d_d_theta_square_length1_dalpha_theta = phi^2*sin(2*alpha);
 d_d_theta_square_length1_dphi_theta = 2*phi*sin(alpha)^2;
 % finallay derivative of 'd_theta_square_length1'
 % with respect to 'vartheta'
d_d_theta_square_length1_vartheta = d_d_theta_square_length1_alpha*...
    dalpha_vartheta + d_d_theta_square_length1_phi*dphi_vartheta +...
    d_d_theta_square_length1_dalpha_theta*d_dalpha_theta_vartheta +...
    d_d_theta_square_length1_dphi_theta*d_phi_theta_vartheta;
% d_theta_square_length2 = 2*...
%     trapz(linspace(0, 1, T), sum((parq0_theta).*(parq0 - q)));
% parallel transport of 'q0' from 'x0' to 'x'
% parq0 = parallel_transport_PL(q0, x0, x, theta);
parq0_vartheta = d_parallel_transport_PL_xu(q0, x0, x, theta, x_vartheta) +...
    d_parallel_transport_PL_theta(q0, x0, x, theta)*theta_vartheta;
% derivative of 'parq0' with respect to 'theta'
% parq0_theta = d_parallel_transport_PL_theta(q0, x0, x, theta);
parq0_theta_vartheta = d2_parallel_transport_PL_xutheta(q0, x0, x, theta, x_vartheta) +...
    d2_parallel_transport_PL_theta(q0, x0, x, theta)*theta_vartheta;
d_d_theta_square_length2_vartheta = 2*...
    trapz(linspace(0, 1, T), sum((parq0_theta_vartheta).*(parq0 - q))) +...
    trapz(linspace(0, 1, T), sum((parq0_theta).*(parq0_vartheta - q_vartheta)));
gradient_vartheta = 2*d_theta_square_length*...
    (d_d_theta_square_length1_vartheta + d_d_theta_square_length2_vartheta);
end

