function [gradient_x, square_length] = DxSquareLength(pi, p, theta)
%%% the derivative 'gradient_x' of the square length 'SquareLength'
%%% with respect to 'x' on the space of 
%%% parameterized curves
%%% Inputs: two curves 'p1' and 'p2', 'theta' which determines 
%%%         the baseline 'beta'
%%% Outputs: the derivative 'gradient_x' 
%%%          of the square length 'SquareLength'
%%%          with respect to 'x'
%%%          and the square length 'square_length'
T = size(pi, 2); % T=size(p2, 2)
% we have starting points on S2
xi = pi(:, 1);
x = p(:, 1);
% and SRVC at tangent spaces
qi = path_to_q(pi);
q = path_to_q(p);
% alpha and phi
alpha = acos(sin(theta)*sqrt((1+dot(xi, x))/2));
phi = 2*asin(sqrt((1-dot(xi, x))/(1+cos(theta)^2 -...
             dot(xi, x)*sin(theta)^2)));
% gradient of alpha and phi with respect to x         
Xix = dot(xi, x);
dalpha_x = -sin(theta)*(xi-dot(xi, x)*x)/...
    (4*(Xix/2 + 1/2)^(1/2)*(1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2)); 
dphi_x = -(2*sin(theta)^2 - 2)*(xi-dot(xi, x)*x)/...
    ((1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
% gradient of length of beta with respect to x
gradient1 = 2*phi*sin(alpha)*...
    (dphi_x*sin(alpha) + phi*dalpha_x*cos(alpha));
% gradient of integral with respect to x
idd = [1 0 0; 0 1 0; 0 0 1];
gradient2 = 0;
parqi = parallel_transport_PL(qi, xi, x, theta);
parqi_x = d_parallel_transport_PL_x(qi, xi, x, theta);
%gradient2 = -trapz(linspace(0,1,T), sum( (parqjn-parqj).*(parqjn-parqj) )  );
for j=1:(T-1)
 gradient2 = gradient2+...
   ((idd-x*x')*parqi_x(:, :, j)'*(parqi(:,j)-q(:, j))+...
    (idd-x*x')*parqi_x(:, :, j+1)'*(parqi(:,j+1)-q(:, j+1)))/((T-1));
end
gradient_x = gradient1 + gradient2;
square_length = SquareLength(pi, p, theta);
end