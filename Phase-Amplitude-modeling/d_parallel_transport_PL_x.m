function dparqi_x = d_parallel_transport_PL_x(qi, xi, x, theta)
% differential of parallel transport of curve 'qi' 
% from the tangent space at 'xi'
% to the tangent space at 'x'
% along the p-optimal curve defined by the angle 'theta'
% with respect to 'x'
T = size(qi, 2);
idd = [1 0 0; 0 1 0; 0 0 1];
phi = 2*asin(sqrt((1-dot(xi, x))/(1+cos(theta)^2-...
             dot(xi, x)*sin(theta)^2)));  
if (phi <= pi/2)
 dparqi_x = d_parallel_transport_PL_x2(qi, xi, x, theta);
else
 normal = (xi + x)*sin(theta)/norm(xi + x)+...
    cross(xi, x)*cos(theta)/norm(cross(xi, x));
 % gradient of 'alpha', 'phi', and 'normal'
 % with respect to x
 Xix = dot(xi, x);
 dphi_x = -(2*sin(theta)^2 - 2)*(xi)/...
    ((1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
 dnormal_x = (sin(theta)/norm(xi+x))*idd -...
   (sin(theta)/norm(xi+x)^3)*(xi+x)*(xi+x)' +...
   (cos(theta)/norm(cross(xi, x)))*ccross(xi) -...
   (cos(theta)/norm(cross(xi, x))^3)*cross(xi, x)*...
   cross(xi, x)'*ccross(xi);
 % middle point on the segment (x1, x2)
 s = 1/2;
 x12 = xi*cos(s*phi) + cross(normal, xi)*sin(s*phi) +...
    normal*dot(normal, xi)*(1 - cos(s*phi));
 % the derivative of the middle point on 
 % the segment (x1, x2) with respect to 'x'
 x12_x = -sin(s*phi)*s*xi*dphi_x' +...
     sin(s*phi)*ccross(-xi)*dnormal_x +...
     cos(s*phi)*s*cross(normal, xi)*dphi_x' +...
     dot(normal, xi)*(1 - cos(s*phi))*dnormal_x +...
     (1 - cos(s*phi))*normal*xi'*dnormal_x +...
     dot(normal, xi)*sin(s*phi)*s*normal*dphi_x';
 % parallel transport of q1 from x1 to x12 along the beta
 xx12 = (xi + x12)/norm(xi + x12);
 % the derivative of 'xx12' with respect to 'x'
 xx12_x = (x12_x)/norm(xi + x12) -...
     (xi + x12)*(xi + x12)'*x12_x/(norm(xi+x12)^3);
 theta11 = asin(dot(xx12, normal));
 % the derivative of 'theta11' with respect to 'x'
 theta11_x = ((normal'*xx12_x) +...
     (xx12'*dnormal_x))'/(1 - dot(xx12, normal)^2)^(1/2);
 parq11 = parallel_transport_PL2(qi, xi, x12, theta11);
 % the derivative of parq11 with respect to theta
 parq11_x2 = d_parallel_transport_PL_x2(qi, xi, x12, theta11);
 parq11_theta = d_parallel_transport_PL_theta2(qi, xi, x12, theta11);
 for j=1:T
  parq11_x2x(:, :, j) = (idd-x12*x12')*parq11_x2(:, :, j)*x12_x;
  parq11_thetax(:, :, j) = (idd-x12*x12')*parq11_theta(:, j)*theta11_x';
 end
 parq11_x = parq11_x2x + parq11_thetax;
 % parallel transport of parq11 from x12 to x2 along the beta 
 xx12 = (x12 + x)/norm(x12 + x);
 % the derivative of xx12 with respect to 'x'
 xx12_x = (x12_x)/norm(x12 + x) -...
     (x12 + x)*(x12 + x)'*x12_x/(norm(x12 + x)^3);
 theta22 = asin(dot(xx12, normal));
 % the derivative of 'theta11' with respect to 'x'
 theta22_x = (normal'*xx12_x +...
     xx12'*dnormal_x)'/(1 - dot(xx12, normal)^2)^(1/2);
 % parqi = parallel_transport_PL2(parq11, x12, x, theta22);
 % the derivative of 'parq11' with respect to 'x'
 parq1_x1 = d_parallel_transport_PL_x1(parq11, x12, x, theta22);
 parq1_theta = d_parallel_transport_PL_theta2(parq11, x12, x, theta22);
 for j=1:T
  parq1_x1x(:, :, j) = (idd-x*x')*parq1_x1(:, :, j)*x12_x;
  parq1_thetax(:, :, j) = (idd-x*x')*parq1_theta(:, j)*theta22_x';
 end
 dparqi_x = parallel_transport_PL22(parq11_x, x12, x, theta22) +...
     d_parallel_transport_PL_x2(parq11, x12, x, theta22) +...
     parq1_x1x + parq1_thetax;
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

function parq111 = parallel_transport_PL22(q1, x1, x2, theta)
% parallel transport of vector 'v1' from the tangent space at 'x1'
% to the tangent space at 'x2'
alpha = acos(sin(theta)*sqrt((1+dot(x1, x2))/2));
phi = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta)^2-...
             dot(x1, x2)*sin(theta)^2)));
normal = (x1+x2)*sin(theta)/norm(x1+x2)+...
     cross(x1, x2)*cos(theta)/norm(cross(x1,x2));
dbeta0 = -x1*phi*sin(0) + cross(normal, x1)*phi*cos(0) +...
       normal*dot(normal, x1)*phi*sin(0); 
if (dot(cross(cross(normal, x1), normal), x1)>0)
 e1=cross(cross(normal, x1), normal)/norm(cross(normal, x1));
else
 e1=-cross(cross(normal, x1), normal)/norm(cross(normal, x1));
end
if (dot(cross(normal, x1), dbeta0)>0)
 e2=cross(normal, x1)/norm(cross(normal, x1));
else
 e2=-cross(normal, x1)/norm(cross(normal, x1)); 
end
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1= -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
v = v1*e1 + v2*e2 + v3*normal;
T=size(q1,2);
parq111=q1;
for j=1:T
 parq111(:, :, j) = v*e2'*q1(:, :, j) +...
     cross(v, x2)*cross(e2, x1)'*q1(:, :, j);
end
end

function dparq1_theta = d_parallel_transport_PL_theta2(q1, x1, x2, theta)
% differential of parallel transport of vector 'v1' from the tangent space at 'x1'
% to the tangent space at 'x2'
% along the p-optimal curve defined by the angle 'theta'
if (x1==x2)
    dparq1_theta=0*q1;
else    
T=size(q1,2);
alpha = acos(sin(theta)*sqrt((1+dot(x1, x2))/2));
phi = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta)^2-...
             dot(x1, x2)*sin(theta)^2)));
normal = (x1+x2)*sin(theta)/norm(x1+x2)+...
     cross(x1,x2)*cos(theta)/norm(cross(x1,x2));
dbeta0 = -x1*phi*sin(0) + cross(normal, x1)*phi*cos(0) +...
       normal*dot(normal, x1)*phi*sin(0);
% derivative of alpha, phi, and normal with respect to theta
Xix = dot(x1, x2);
dalpha_theta = -(cos(theta)*(Xix/2 + 1/2)^(1/2))/...
    (1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2); 
dphi_theta = -(sin(2*theta)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2))/...
    ((1 - sin(theta)^2)^(1/2)*(sin(theta)^2 + Xix*sin(theta)^2 - 2));
dnormal_theta = (x1+x2)*cos(theta)/norm(x1+x2) -...
     cross(x1,x2)*sin(theta)/norm(cross(x1,x2));
% ee2 and dee2_theta
ee2=cross(normal, x1)/norm(cross(normal, x1));
dee2_theta=cross(dnormal_theta, x1)/norm(cross(normal, x1)) -...
     cross(normal, x1)*...
     dot(cross(dnormal_theta, x1), cross(normal, x1))/...
     norm(cross(normal, x1))^3;
% e2 and de2_theta 
if (dot(cross(normal, x1), dbeta0)>0)
    e2=ee2;
    de2_theta=dee2_theta;
else
    e2=-ee2;
    de2_theta=-dee2_theta;
end
% e1 and de1_theta
if (dot(cross(cross(normal, x1), normal), x1)>0)
    e1=cross(ee2, normal);
    de1_theta=cross(dee2_theta, normal)+cross(ee2, dnormal_theta);
else
    e1=-cross(ee2, normal);
    de1_theta=-cross(dee2_theta, normal)-cross(ee2, dnormal_theta);
end    
% v1, v2, and v3
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% derivative v1, v2, and v3 with respect to theta
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
dv3_theta=dv3_alpha*dalpha_theta+dv3_phi*dphi_theta;
dv2_theta=dv2_alpha*dalpha_theta+dv2_phi*dphi_theta;
dv1_theta=dv1_alpha*dalpha_theta+dv1_phi*dphi_theta;
v = v1*e1 + v2*e2 + v3*normal;
dv_theta = dv1_theta*e1 + dv2_theta*e2 + dv3_theta*normal +...
            +v1*de1_theta + v2*de2_theta + v3*dnormal_theta;
dparq1_theta=q1;
for j=1:T
   %{ 
     parq1(:, j) = dot(e2, q1(:, j))*vm +...
     dot(cross(e2, xi), q1(:, j))*cross(vm, x);
   %}
 dparq1_theta(:, j) = dot(e2, q1(:, j))*dv_theta+...
     dot(cross(e2, x1), q1(:, j))*cross(dv_theta, x2)+...
     dot(de2_theta, q1(:, j))*v+...
     dot(cross(de2_theta, x1), q1(:, j))*cross(v, x2);
end
end
end

function dparqi_xi = d_parallel_transport_PL_x1(qi, xi, x, theta)
% differential of parallel transport of curve 'qi' 
% from the tangent space at 'xi'
% to the tangent space at 'x'
% along the p-optimal curve defined by the angle 'theta'
% with respect to 'x'
T = size(qi,2);
% alpha, phi, normal, and ee2
alpha = acos(sin(theta)*sqrt((1+dot(xi, x))/2));
phi = 2*asin(sqrt((1-dot(xi,x))/(1+cos(theta)^2-...
             dot(xi, x)*sin(theta)^2)));
normal = (xi+x)*sin(theta)/norm(xi+x)+...
     cross(xi,x)*cos(theta)/norm(cross(xi, x));
ee2 = cross(normal, xi)/norm(cross(normal, xi));

dbeta0 = -xi*phi*sin(0) + cross(normal, xi)*phi*cos(0) +...
       normal*dot(normal, xi)*phi*sin(0);
% gradient of alpha, phi, normal, and e2
Xix = dot(xi, x);
dalpha_xi = -sin(theta)*(x)/...
    (4*(Xix/2 + 1/2)^(1/2)*(1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2)); 
dphi_xi = -(2*sin(theta)^2 - 2)*(x)/...
    ((1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
idd = [1 0 0; 0 1 0; 0 0 1];
dnormal_xi = (sin(theta)/norm(xi+x))*idd -...
   (sin(theta)/norm(xi+x)^3)*(xi+x)*(xi+x)' -...
   (cos(theta)/norm(cross(xi, x)))*ccross(x) +...
   (cos(theta)/norm(cross(xi, x))^3)*cross(xi, x)*...
   cross(xi, x)'*ccross(x);
dee2_xi = cross(normal, xi)*cross(normal, xi)'*...
        (ccross(xi)*dnormal_xi - ccross(normal))/...
        norm(cross(normal, xi))^3 +...
        (ccross(normal)-ccross(xi)*dnormal_xi)/norm(cross(normal, xi));
% differential of e2 and e3 with respect to 'x'
if (dot(cross(normal, xi), dbeta0)>0)
   e2 = ee2;
   de2_xi = dee2_xi;
else
   e2 = -ee2;
   de2_xi = -dee2_xi;
end
if (dot(cross(cross(normal, xi), normal), xi)>0)
    e1 = cross(ee2, normal);
    de1_xi = -ccross(normal)*dee2_xi + ccross(ee2)*dnormal_xi;
else
    e1 = -cross(ee2, normal);
    de1_xi = ccross(normal)*dee2_xi - ccross(ee2)*dnormal_xi;
end
% v1, v2, and v3
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% gradient v1, v2, and v3
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
dv3_xi = dv3_alpha*dalpha_xi + dv3_phi*dphi_xi;
dv2_xi = dv2_alpha*dalpha_xi + dv2_phi*dphi_xi;
dv1_xi = dv1_alpha*dalpha_xi + dv1_phi*dphi_xi;
% v
v = v1*e1 + v2*e2 + v3*normal;
% differential of v with respect to 'x'
dv_xi = e1*dv1_xi' + e2*dv2_xi' + normal*dv3_xi' +...
            +v1*de1_xi + v2*de2_xi + v3*dnormal_xi;
% differential of parallel transport of qi with respect to x        
dparqi_xi = zeros(3, 3, T);
for j=1:T
   %{ 
     parqi(:, j) = dot(e2, qi(:, j))*v +...
     dot(cross(e2, xi), qi(:, j))*cross(v, x);
   %}
 dparqi_xi(:, :, j) = v*qi(:, j)'*de2_xi -...
     cross(v, x)*qi(:, j)'*ccross(xi)*de2_xi +...
     cross(v, x)*qi(:, j)'*ccross(e2) +...
     dot(e2, qi(:, j))*dv_xi -...
     dot(cross(e2, xi), qi(:, j))*ccross(x)*dv_xi;
end
end

function dparqi_x = d_parallel_transport_PL_x2(qi, xi, x, theta)
% differential of parallel transport of curve 'qi' 
% from the tangent space at 'xi'
% to the tangent space at 'x'
% along the p-optimal curve defined by the angle 'theta'
% with respect to 'x'
T=size(qi,2);
% alpha, phi, normal, and ee2
alpha = acos(sin(theta)*sqrt((1+dot(xi, x))/2));
phi = 2*asin(sqrt((1-dot(xi,x))/(1+cos(theta)^2-...
             dot(xi, x)*sin(theta)^2)));
normal = (xi+x)*sin(theta)/norm(xi+x)+...
     cross(xi,x)*cos(theta)/norm(cross(xi, x));
ee2 = cross(normal, xi)/norm(cross(normal, xi));

dbeta0 = -xi*phi*sin(0) + cross(normal, xi)*phi*cos(0) +...
       normal*dot(normal, xi)*phi*sin(0);
% gradient of alpha, phi, normal, and e2
Xix = dot(xi, x);
dalpha_x = -sin(theta)*(xi)/...
    (4*(Xix/2 + 1/2)^(1/2)*(1 - sin(theta)^2*(Xix/2 + 1/2))^(1/2)); 
dphi_x = -(2*sin(theta)^2 - 2)*(xi)/...
    ((1 - sin(theta)^2)^(1/2)*(1 - Xix)^(1/2)*(Xix + 1)^(1/2)*...
    (sin(theta)^2 + Xix*sin(theta)^2 - 2));
idd = [1 0 0; 0 1 0; 0 0 1];
dnormal_x = (sin(theta)/norm(xi+x))*idd -...
   (sin(theta)/norm(xi+x)^3)*(xi+x)*(xi+x)' +...
   (cos(theta)/norm(cross(xi, x)))*ccross(xi) -...
   (cos(theta)/norm(cross(xi, x))^3)*cross(xi, x)*...
   cross(xi, x)'*ccross(xi);
dee2_x = cross(normal, xi)*cross(normal, xi)'*ccross(xi)*dnormal_x/...
        norm(cross(normal, xi))^3-...
        ccross(xi)*dnormal_x/norm(cross(normal, xi));
% differential of e2 and e3 with respect to 'x'
if (dot(cross(normal, xi), dbeta0)>=0)
   e2 = ee2;
   de2_x = dee2_x;
else
   e2 = -ee2;
   de2_x = -dee2_x;
end
if (dot(cross(cross(normal, xi), normal), xi)>=0)
    e1 = cross(ee2, normal);
    de1_x = -ccross(normal)*dee2_x + ccross(ee2)*dnormal_x;
else
    e1 = -cross(ee2, normal);
    de1_x = ccross(normal)*dee2_x - ccross(ee2)*dnormal_x;
end
% v1, v2, and v3
v3 = -sin(alpha)*sin(phi*cos(alpha));
v2 = cos(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*sin(phi)*sin(phi*cos(alpha));
v1 = -sin(phi)*cos(phi*cos(alpha))+...
    cos(alpha)*cos(phi)*sin(phi*cos(alpha));
% gradient v1, v2, and v3
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
dv3_x = dv3_alpha*dalpha_x + dv3_phi*dphi_x;
dv2_x = dv2_alpha*dalpha_x + dv2_phi*dphi_x;
dv1_x = dv1_alpha*dalpha_x + dv1_phi*dphi_x;
% v
v = v1*e1 + v2*e2 + v3*normal;
% differential of v with respect to 'x'
dv_x = e1*dv1_x' + e2*dv2_x' + normal*dv3_x' +...
            +v1*de1_x + v2*de2_x + v3*dnormal_x;
% differential of parallel transport of qi with respect to x        
dparqi_x = zeros(3, 3, T);
for j=1:T
   %{ 
     parqi(:, j) = dot(e2, qi(:, j))*v +...
     dot(cross(e2, xi), qi(:, j))*cross(v, x);
   %}
 dparqi_x(:, :, j) = v*qi(:, j)'*de2_x -...
     cross(v, x)*qi(:, j)'*ccross(xi)*de2_x +...
     dot(e2, qi(:, j))*dv_x -...
     dot(cross(e2, xi), qi(:, j))*ccross(x)*dv_x +...
     dot(cross(e2, xi), qi(:, j))*ccross(v);
end
end

function ma = ccross(a)
ma=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end