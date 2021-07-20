function visparq1 = visparallel_transport_PL(q1, x1, x2, N2, theta)
% parallel transport of vector 'v1' from the tangent space at 'x1'
% to the tangent space at 'x2'
% along the p-optimal curve defined by the angle 'theta'
alpha = acos(sin(theta)*sqrt((1+dot(x1, x2))/2));
phi = 2*asin(sqrt((1-dot(x1,x2))/(1+cos(theta)^2-...
             dot(x1, x2)*sin(theta)^2)));         
normal = (x1+x2)*sin(theta)/norm(x1+x2)+...
     cross(x1,x2)*cos(theta)/norm(cross(x1,x2));
%T=size(q1, 2);
s2=linspace(0,1,N2);
for j=1:N2
beta2(:, j) = cross(cross(normal, x1), normal)*cos(s2(j)*phi)*sin(alpha)/...
     norm(cross(cross(normal, x1),normal)) +...
     cross(normal, x1)*sin(s2(j)*phi)*sin(alpha)/norm(cross(normal, x1)) +...
     normal*cos(alpha);
end
for j=1:N2
  xx12 = (x1 + beta2(:, j))/norm(x1 + beta2(:, j));
  theta1 = asin(dot(xx12, normal));
  visparq1(:, :, j) = parallel_transport_PL(q1, x1, beta2(:, j), theta1);
end
end