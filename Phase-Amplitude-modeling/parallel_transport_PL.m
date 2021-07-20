function parq1 = parallel_transport_PL(q1, x1, x2, theta)
% parallel transport of vector 'v1' from the tangent space at 'x1'
% to the tangent space at 'x2'
 if (norm(x1 - x2) < 0.00000000001)
  parq1 = q1;
 else 
  phi = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta)^2-...
             dot(x1, x2)*sin(theta)^2)));  
  if (phi < pi/2)
    parq1 = parallel_transport_PL2(q1, x1, x2, theta);
  else
    normal = (x1 + x2)*sin(theta)/norm(x1 + x2) +...
          cross(x1, x2)*cos(theta)/norm(cross(x1, x2));
  % middle point on the segment (x1, x2)
    s = 1/2;
    x12 = x1*cos(s*phi) + cross(normal, x1)*sin(s*phi) +...
          normal*dot(normal, x1)*(1 - cos(s*phi));
  % parallel transport of q1 from x1 to x12 along the beta
    xx12 = (x1 + x12)/norm(x1 + x12);
    theta1 = asin(dot(xx12, normal));
    parq11 = parallel_transport_PL2(q1, x1, x12, theta1);
  % parallel transport of parq11 from x12 to x2 along the beta 
    xx12 = (x12 + x2)/norm(x12 + x2);
    theta2 = asin(dot(xx12, normal));
    parq1 = parallel_transport_PL2(parq11, x12, x2, theta2);
  end
 end
end
function parq111 = parallel_transport_PL2(q1, x1, x2, theta)
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
 parq111(:, j) = dot(e2, q1(:, j))*v +...
     dot(cross(e2, x1), q1(:, j))*cross(v, x2);
end
end


