function dparq1_theta = d_parallel_transport_PL_theta(q1, x1, x2, theta)
phi = 2*asin(sqrt((1 - dot(x1, x2))/(1 + cos(theta)^2 -...
             dot(x1, x2)*sin(theta)^2)));  
if (phi <= pi/2)
 dparq1_theta = d_parallel_transport_PL_theta2(q1, x1, x2, theta);
else
 normal = (x1 + x2)*sin(theta)/norm(x1 + x2) +...
          cross(x1, x2)*cos(theta)/norm(cross(x1, x2));
 % derivative of 'phi' and 'normal' with respect to 'theta'
 Xix = dot(x1, x2); 
 dphi_theta = -(sin(2*theta)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2))/...
    ((1 - sin(theta)^2)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
 dnormal_theta = (x1 + x2)*cos(theta)/norm(x1 + x2) -...
     cross(x1, x2)*sin(theta)/norm(cross(x1, x2));
 % middle point on the segment (x1, x2)
 s = 1/2;
 x12 = x1*cos(s*phi) + cross(normal, x1)*sin(s*phi) +...
    normal*dot(normal, x1)*(1 - cos(s*phi));
 % the derivative of the middle point on the segment (x1, x2) 
% with respect to 'theta'
 x12_theta = -x1*sin(s*phi)*s*dphi_theta +...% the second term
     cross(dnormal_theta, x1)*sin(s*phi) +...
     cross(normal, x1)*cos(s*phi)*s*dphi_theta +...% the third term
     dnormal_theta*dot(normal, x1)*(1 - cos(s*phi)) +...
     normal*dot(dnormal_theta, x1)*(1 - cos(s*phi)) +...
     normal*dot(normal, x1)*sin(s*phi)*s*dphi_theta;
 % parallel transport of 'qi' from 'x1' to 'x12' along 
 % the 'beta11' which is the first half of 'beta'
 % defined by angle between 'normal' and 'xx12'
 xx12 = (x1 + x12)/norm(x1 + x12);
 % the derivative of 'xx12' with respect to 'theta'
 xx12_theta = (x12_theta)/norm(x1 + x12) -...
     (x1 + x12)*dot(x1 + x12, x12_theta)/(norm(x1 + x12)^3);
 % the anlge 'theta11' defining the 'beta11'
 theta11 = asin(dot(xx12, normal));
 % the derivative of 'theta11' with respect to 'theta'
 theta11_theta = (dot(xx12_theta, normal) +...
     dot(xx12, dnormal_theta))/(1 - dot(xx12, normal)^2)^(1/2);
 % the parallel tranport of 'q1' along the 'beta11'
 parq11 = parallel_transport_PL2(q1, x1, x12, theta11);
 % the derivative of 'parq11' with respect to 'theta'
 parq11_theta = d_parallel_transport_PL_x222(q1, x1, x12, theta11, x12_theta) +...
     d_parallel_transport_PL_theta2(q1, x1, x12, theta11)*theta11_theta;
 % parallel transport of parq11 from x12 to x2 along the beta 
 xx12 = (x12 + x2)/norm(x12 + x2);
 % the derivative of 'xx12' with respect to 'theta'
 xx12_theta = (x12_theta)/norm(x12 + x2) -...
     (x12 + x2)*dot(x12 + x2, x12_theta)/(norm(x12 + x2)^3);
 % the anlge 'theta22' defining the 'beta22'
 theta22 = asin(dot(xx12, normal));
 % the derivative of 'theta22' with respect to 'theta'
 theta22_theta = (dot(xx12_theta, normal) +...
     dot(xx12, dnormal_theta))/(1 - dot(xx12, normal)^2)^(1/2);
 % the final parallel transport of 'parq11' from 'x12' to 'x2'
 % parq1 = parallel_transport_PL2(parq11, x12, x2, theta22);
 % the derivative of parq11 with respect to 'theta'
 dparq1_theta = parallel_transport_PL2(parq11_theta, x12, x2, theta22) +...
     d_parallel_transport_PL_x122(parq11, x12, x2, theta22, x12_theta) +...
     d_parallel_transport_PL_theta2(parq11, x12, x2, theta22)*theta22_theta;
end
end

function parq111 = parallel_transport_PL2(q1, x1, x2, theta)
% parallel transport of vector 'q1'
% from the tangent space at 'x1'
% to the tangent space at 'x2'
% along the p-optimal baseline 'beta'
% determined by the angle 'theta'
alpha = acos(sin(theta)*sqrt((1 + dot(x1, x2))/2));
phi = 2*asin(sqrt((1 - dot(x1, x2))/(1 + cos(theta)^2 -...
             dot(x1, x2)*sin(theta)^2)));
normal = (x1 + x2)*sin(theta)/norm(x1 + x2) +...
     cross(x1, x2)*cos(theta)/norm(cross(x1, x2));
% direction of the p-optimal curve at the beginning
dbeta0 = -x1*phi*sin(0) + cross(normal, x1)*phi*cos(0) +...
       normal*dot(normal, x1)*phi*sin(0);
if (dot(cross(cross(normal, x1), normal), x1)>0)
 e1 = cross(cross(normal, x1), normal)/norm(cross(normal, x1));
else
 e1 = -cross(cross(normal, x1), normal)/norm(cross(normal, x1));
end
if (dot(cross(normal, x1), dbeta0)>0)
 e2 = cross(normal, x1)/norm(cross(normal, x1));
else
 e2 = -cross(normal, x1)/norm(cross(normal, x1)); 
end
% coefficient of vector 'v' (parallel transport of 'e2')
% with basis 'e1', 'e2', and 'normal'
% 'v1', 'v2', and 'v3'
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1= -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% the vector 'v'
v = v1*e1 + v2*e2 + v3*normal;
T = size(q1, 2);
parq111 = q1;
for j=1:T
 parq111(:, j) = dot(e2, q1(:, j))*v +...
     dot(cross(e2, x1), q1(:, j))*cross(v, x2);
end
end

% the first derivatives

function dparq1_theta = d_parallel_transport_PL_theta2(q1, x1, x2, theta)
% differential of parallel transport of vector 'v1' from the tangent space at 'x1'
% to the tangent space at 'x2'
% along the p-optimal curve defined by the angle 'theta'
if (x1 == x2)
    dparq1_theta = 0*q1;
else    
T = size(q1, 2);
% 'alpha', 'phi', and 'normal'
alpha = acos(sin(theta)*sqrt((1+dot(x1, x2))/2));
phi = 2*asin(sqrt((1 - dot(x1, x2))/(1 + cos(theta)^2 -...
             dot(x1, x2)*sin(theta)^2)));
normal = (x1 + x2)*sin(theta)/norm(x1 + x2) +...
     cross(x1, x2)*cos(theta)/norm(cross(x1, x2));
% direction of the p-optimal curve at the beginning
dbeta0 = -x1*phi*sin(0) + cross(normal, x1)*phi*cos(0) +...
       normal*dot(normal, x1)*phi*sin(0);
% derivatives of 'alpha', 'phi', and 'normal' with respect to 'theta'
Xix = dot(x1, x2);
dalpha_theta = -(cos(theta)*(Xix/2 + 1/2)^(1/2))/...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2); 
dphi_theta = -(sin(2*theta)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2))/...
    ((1 - sin(theta)^2)^(1/2)*(sin(theta)^2 + Xix*sin(theta)^2 - 2));
dnormal_theta = (x1+x2)*cos(theta)/norm(x1+x2) -...
     cross(x1, x2)*sin(theta)/norm(cross(x1,x2));
% 'ee2' 
ee2 = cross(normal, x1)/norm(cross(normal, x1));
% the derivative of 'ee2' with respect to 'theta'
dee2_theta = cross(dnormal_theta, x1)/norm(cross(normal, x1)) -...
     cross(normal, x1)*...
     dot(cross(dnormal_theta, x1), cross(normal, x1))/...
     norm(cross(normal, x1))^3;
% 'e2' and the derivative of 'e2' with respect to 'theta' 
if (dot(cross(normal, x1), dbeta0) > 0)
    e2 = ee2;
    de2_theta = dee2_theta;
else
    e2 = -ee2;
    de2_theta = -dee2_theta;
end
% 'e1' and the derivative of 'e1' with respect to 'theta' 
if (dot(cross(cross(normal, x1), normal), x1)>0)
    e1 = cross(ee2, normal);
    de1_theta = cross(dee2_theta, normal) +...
        cross(ee2, dnormal_theta);
else
    e1 = -cross(ee2, normal);
    de1_theta = -cross(dee2_theta, normal) -...
        cross(ee2, dnormal_theta);
end    
% coefficient of vector 'v' (parallel transport of 'e2')
% with basis 'e1', 'e2', and 'normal'
% 'v1', 'v2', and 'v3'
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% derivatives of 'v1', 'v2', and 'v3' with respect to
% 'alpha' and 'phi'
dv3_alpha = phi*sin(alpha)^2*cos(phi*cos(alpha)) -...
    sin(phi*cos(alpha))*cos(alpha);
dv3_phi = -cos(alpha)*sin(alpha)*cos(phi*cos(alpha));
dv2_alpha = -sin(alpha)*(sin(phi*cos(alpha))*sin(phi) -...
    phi*sin(phi*cos(alpha))*cos(phi) +...
    phi*cos(alpha)*sin(phi)*cos(phi*cos(alpha)));
dv2_phi = sin(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
dv1_alpha = -sin(alpha)*(sin(phi*cos(alpha))*cos(phi) +...
    phi*sin(phi*cos(alpha))*sin(phi) +...
    phi*cos(alpha)*cos(phi)*cos(phi*cos(alpha)));
dv1_phi = cos(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
% derivatives of 'v1', 'v2', and 'v3' with respect to 'theta'
dv3_theta = dv3_alpha*dalpha_theta + dv3_phi*dphi_theta;
dv2_theta = dv2_alpha*dalpha_theta + dv2_phi*dphi_theta;
dv1_theta = dv1_alpha*dalpha_theta + dv1_phi*dphi_theta;
% the vector 'v'
v = v1*e1 + v2*e2 + v3*normal;
% the derivative of 'v' with respect to 'theta'
dv_theta = dv1_theta*e1 + dv2_theta*e2 + dv3_theta*normal +...
            v1*de1_theta + v2*de2_theta + v3*dnormal_theta;
dparq1_theta = q1;
for j=1:T
   %{ 
     parq1(:, j) = dot(e2, q1(:, j))*v +...
     dot(cross(e2, x1), q1(:, j))*cross(v, x2);
   %}
 dparq1_theta(:, j) = dot(de2_theta, q1(:, j))*v +...
     dot(e2, q1(:, j))*dv_theta +...
     dot(cross(de2_theta, x1), q1(:, j))*cross(v, x2) +...
     dot(cross(e2, x1), q1(:, j))*cross(dv_theta, x2);
end
end
end

function dparqi_xiu = d_parallel_transport_PL_x122(qi, xi, x, theta, uu)
% differential of parallel transport of curve 'qi' 
% from the tangent space at 'xi'
% to the tangent space at 'x'
% along the p-optimal curve 'beta'
% defined by the angle 'theta'
% with respect to 'xi' 
% along the direction 'uu'
T = size(qi, 2);
% 'alpha', 'phi', 'normal', and 'ee2'
alpha = acos(sin(theta)*sqrt((1 + dot(xi, x))/2));
phi = 2*asin(sqrt((1 - dot(xi, x))/(1 + cos(theta)^2 -...
             dot(xi, x)*sin(theta)^2)));
normal = (xi + x)*sin(theta)/norm(xi + x) +...
     cross(xi, x)*cos(theta)/norm(cross(xi, x));
ee2 = cross(normal, xi)/norm(cross(normal, xi));
% direction of the p-optimal curve 'beta' at the beginning
dbeta0 = -xi*phi*sin(0) + cross(normal, xi)*phi*cos(0) +...
       normal*dot(normal, xi)*phi*sin(0);
% differentials of 'alpha', 'phi', 'normal', and 'ee2'
% with respect to 'xi' along the direction 'uu'
Xix = dot(xi, x);
dalpha_xiu = -sin(theta)*dot(uu, x)/...
    (4*(Xix/2 + 1/2)^(1/2)*(1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2)); 
dphi_xiu = -(2*sin(theta)^2 - 2)*dot(uu, x)/...
    ((1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
dnormal_xiu = uu*(sin(theta)/norm(xi + x)) -...
   (xi + x)*(sin(theta)/norm(xi + x)^3)*dot(xi + x, uu) +...% the second term
   cross(uu, x)*(cos(theta)/norm(cross(xi, x))) -...
   cross(xi, x)*(cos(theta)/norm(cross(xi, x))^3)*...
   dot(cross(xi, x), cross(uu, x));
dee2_xiu = (cross(dnormal_xiu, xi) + cross(normal, uu))/...
    norm(cross(normal, xi)) -...
    cross(normal, xi)*...
    dot(cross(normal, xi), cross(dnormal_xiu, xi) +...
    cross(normal, uu))/norm(cross(normal, xi))^3;
% differentials of 'e2' and 'e1' with respect to 'xi'
% along the direction 'uu'
if (dot(cross(normal, xi), dbeta0) > 0)
   e2 = ee2;
   de2_xiu = dee2_xiu;
else
   e2 = -ee2;
   de2_xiu = -dee2_xiu;
end
if (dot(cross(cross(normal, xi), normal), xi) > 0)
    e1 = cross(ee2, normal);
    de1_xiu = cross(dee2_xiu, normal) + cross(ee2, dnormal_xiu);
else
    e1 = -cross(ee2, normal);
    de1_xiu = -cross(dee2_xiu, normal) - cross(ee2, dnormal_xiu);
end
% coefficient of vector 'v' (parallel transport of 'e2')
% with basis 'e1', 'e2', and 'normal'
% 'v1', 'v2', and 'v3'
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% differentials of 'v1', 'v2', and 'v3' with respect to
% 'alpha' and 'phi'
dv3_alpha = phi*sin(alpha)^2*cos(phi*cos(alpha)) -...
    sin(phi*cos(alpha))*cos(alpha);
dv3_phi = -cos(alpha)*sin(alpha)*cos(phi*cos(alpha));
dv2_alpha = -sin(alpha)*(sin(phi*cos(alpha))*sin(phi) -...
    phi*sin(phi*cos(alpha))*cos(phi) +...
    phi*cos(alpha)*sin(phi)*cos(phi*cos(alpha)));
dv2_phi = sin(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
dv1_alpha = -sin(alpha)*(sin(phi*cos(alpha))*cos(phi) +...
    phi*sin(phi*cos(alpha))*sin(phi) +...
    phi*cos(alpha)*cos(phi)*cos(phi*cos(alpha)));
dv1_phi = cos(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
% differentials of 'v1', 'v2', and 'v3' with respect to 'xi'
% along the direction 'uu'
dv3_xiu = dv3_alpha*dalpha_xiu + dv3_phi*dphi_xiu;
dv2_xiu = dv2_alpha*dalpha_xiu + dv2_phi*dphi_xiu;
dv1_xiu = dv1_alpha*dalpha_xiu + dv1_phi*dphi_xiu;
% the vector 'v'
v = v1*e1 + v2*e2 + v3*normal;
% differential of 'v' with respect to 'xi'
% along the direction 'uu'
dv_xiu = dv1_xiu*e1 + dv2_xiu*e2 + dv3_xiu*normal +...
            +v1*de1_xiu + v2*de2_xiu + v3*dnormal_xiu;
% differential of parallel transport of 'qi' with respect to 'xi'
% along the direction 'uu'
dparqi_xiu = zeros(3, T);
for j=1:T
   %{ 
     parqi(:, j) = dot(e2, qi(:, j))*v +...
     dot(cross(e2, xi), qi(:, j))*cross(v, x);
   %}
 dparqi_xiu(:, j) = dot(de2_xiu, qi(:, j))*v +...
     dot(e2, qi(:, j))*dv_xiu +...% the second term
     dot(cross(de2_xiu, xi), qi(:, j))*cross(v, x) +...
     dot(cross(e2, uu), qi(:, j))*cross(v, x) +...
     dot(cross(e2, xi), qi(:, j))*cross(dv_xiu, x);
end
end

function dparqi_xu = d_parallel_transport_PL_x222(qi, xi, x, theta, uu)
% differential of parallel transport of curve 'qi' 
% from the tangent space at 'xi'
% to the tangent space at 'x'
% along the p-optimal curve 
% defined by the angle 'theta'
% with respect to 'x'
% along the direction 'uu'
T = size(qi, 2);
% 'alpha', 'phi', 'normal', and 'ee2'
alpha = acos(sin(theta)*sqrt((1 + dot(xi, x))/2));
phi = 2*asin(sqrt((1 - dot(xi, x))/(1 + cos(theta)^2 -...
             dot(xi, x)*sin(theta)^2)));
normal = (xi + x)*sin(theta)/norm(xi + x) +...
     cross(xi, x)*cos(theta)/norm(cross(xi, x));
ee2 = cross(normal, xi)/norm(cross(normal, xi));
% direction of the p-optimal curve 'beta' at the beginning
dbeta0 = -xi*phi*sin(0) + cross(normal, xi)*phi*cos(0) +...
       normal*dot(normal, xi)*phi*sin(0);
% differentials of 'alpha', 'phi', 'normal', and 'ee2'
% with respect to 'x'
% along the direction 'uu'
Xix = dot(xi, x);
dalpha_xu = -sin(theta)*dot(xi, uu)/...
    (4*(Xix/2 + 1/2)^(1/2)*...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2)); 
dphi_xu = -(2*sin(theta)^2 - 2)*dot(xi, uu)/...
    ((1 - sin(theta)^2)^(1/2)*...
    (1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
dnormal_xu = uu*(sin(theta)/norm(xi + x)) -...
   (xi + x)*(sin(theta)/norm(xi + x)^3)*dot(xi + x, uu) +...% the second term
   cross(xi, uu)*(cos(theta)/norm(cross(xi, x))) -...
   cross(xi, x)*(cos(theta)/norm(cross(xi, x))^3)*...
   dot(cross(xi, x), cross(xi, uu));
dee2_xu = cross(dnormal_xu, xi)/norm(cross(normal, xi)) -...
    cross(normal, xi)*...
    dot(cross(normal, xi), cross(dnormal_xu, xi))/...
    norm(cross(normal, xi))^3;
% differentials of 'e2' and 'e1' with respect to 'x'
% along the direction 'uu'
if (dot(cross(normal, xi), dbeta0) > 0)
   e2 = ee2;
   de2_xu = dee2_xu;
else
   e2 = -ee2;
   de2_xu = -dee2_xu;
end
if (dot(cross(cross(normal, xi), normal), xi) > 0)
    e1 = cross(ee2, normal);
    de1_xu = cross(dee2_xu, normal) + cross(ee2, dnormal_xu);
else
    e1 = -cross(ee2, normal);
    de1_xu = -cross(dee2_xu, normal) - cross(ee2, dnormal_xu);
end
% coefficient of vector 'v' (parallel transport of 'e2')
% with basis 'e1', 'e2', and 'normal'
% 'v1', 'v2', and 'v3'
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% differentials of 'v1', 'v2', and 'v3' with respect to
% 'alpha' and 'phi'
dv3_alpha = phi*sin(alpha)^2*cos(phi*cos(alpha)) -...
    sin(phi*cos(alpha))*cos(alpha);
dv3_phi = -cos(alpha)*sin(alpha)*cos(phi*cos(alpha));
dv2_alpha = -sin(alpha)*(sin(phi*cos(alpha))*sin(phi) -...
    phi*sin(phi*cos(alpha))*cos(phi) +...
    phi*cos(alpha)*sin(phi)*cos(phi*cos(alpha)));
dv2_phi = sin(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
dv1_alpha = -sin(alpha)*(sin(phi*cos(alpha))*cos(phi) +...
    phi*sin(phi*cos(alpha))*sin(phi) +...
    phi*cos(alpha)*cos(phi)*cos(phi*cos(alpha)));
dv1_phi = cos(phi)*cos(phi*cos(alpha))*(cos(alpha)^2 - 1);
% differentials of 'v1', 'v2', and 'v3' with respect to 'x'
% along the direction 'uu'
dv3_xu = dv3_alpha*dalpha_xu + dv3_phi*dphi_xu;
dv2_xu = dv2_alpha*dalpha_xu + dv2_phi*dphi_xu;
dv1_xu = dv1_alpha*dalpha_xu + dv1_phi*dphi_xu;
% the vector 'v'
v = v1*e1 + v2*e2 + v3*normal;
% differential of 'v' with respect to 'x'
% along the direction 'uu'
dv_xu = dv1_xu*e1 + dv2_xu*e2 + dv3_xu*normal +...
            +v1*de1_xu + v2*de2_xu + v3*dnormal_xu;
% differential of parallel transport of 'qi' with respect to 'x'
% along the direction 'uu'        
dparqi_xu = zeros(3, T);
for j=1:T
   %{ 
     parqi(:, j) = dot(e2, qi(:, j))*v +...
     dot(cross(e2, xi), qi(:, j))*cross(v, x);
   %}
 dparqi_xu(:, j) = dot(de2_xu, qi(:, j))*v +...
    dot(e2, qi(:, j))*dv_xu +...% the second term
    dot(cross(de2_xu, xi), qi(:, j))*cross(v, x) +...
    dot(cross(e2, xi), qi(:, j))*cross(dv_xu, x) +...
    dot(cross(e2, xi), qi(:, j))*cross(v, uu);
end
end