function [p, theta, vartheta, gradient_theta2] = exponential_map(p0, u, w, step, epsilon, NN1, NN2, NN3)
%%% the Exponential map on the space of parameterized curves
%%% the function is finding 'vartheta' that determines the baseline which generates
%%% the geodesic between p1 and p2, which is determined by also 'u' and 'w',
%%% on and the space of parameterized curves
%%% the inputs: the curves 'p0', the direction 'u' of the baseline 'beta' 
%%%             the differece 'w' of TSRVCs,
%%%             the degree of accurace 'epsilon', the step size 'step',
%%%             the maximum number of iterations 'NN1',
%%%             the maximum number of iterations when we don't perform
%%%             steps 'NN2' and 'NN3'
%%% the outputs: the 'vartheta' that determines baseline 'beta',
%%%              the 'theta' that is determined by 'vartheta'
%%%              the square gradient 'gradient_theta2' with respect to 'theta',
%%%              the curve 'p'
% the initial value
 %vartheta = pi/2;
 M = 8;
 t = linspace(asin(norm(u)/(pi)) + 1/1000, pi - asin(norm(u)/(pi)) - 1/1000, M);
 for k=1:M
  grad_theta2(k) =   SquareDthetaSquareLength(p0, u, w, t(k));
 end
 [mingrad_theta2 , mink] = min(grad_theta2);
 vartheta = t(mink);
 [gradient_vartheta, gradient_theta2] = DvarthetaSquareDthetaSquareLength(p0, u, w, vartheta);
 N1 = 0;
 N2 = 0;
 N3 = 0;
 step1 = step; nn = 1;
 while(sqrt(gradient_theta2) > epsilon)
  if (N1 <= NN1)
   N1 = N1 + 1;
   vartheta1 = vartheta - step1*gradient_vartheta;
   if ((vartheta1 > asin(norm(u)/(pi))) && (vartheta1 < pi - asin(norm(u)/(pi))))
    N3 = 0;   
    gradient_theta21 = SquareDthetaSquareLength(p0, u, w, vartheta1);
    if (gradient_theta21 < gradient_theta2)
     N2 = 0;   
     vartheta = vartheta1; 
     [gradient_vartheta, gradient_theta2] = DvarthetaSquareDthetaSquareLength(p0, u, w, vartheta);
    else
     if(N2 > NN2)
      nn = norm(step1*gradient_vartheta);
      X = sprintf('we are very close in neighborhood of %e, \n with the gradient %e \n' , nn, sqrt(gradient_theta2));
      disp(X);
      break
     else
      N2 = N2 + 1; nn = norm(step1*gradient_vartheta);
      step1 = step1/2;
     end
    end
   else
    if (N3 > NN3)
     X = sprintf('we reached the bound 0 or \pi, \n with the gradient %e \n', sqrt(gradient_theta2));
     disp(X);
     break
    else
     N3 = N3 + 1;
     step1 = step1/2;
    end
   end 
  else
   X = sprintf('we reached the number of iterations %d, with the gradient %e \n and the neighborhood %e \n', ...
       NN1, sqrt(gradient_theta2), nn);
   disp(X);
   break
  end
 end 
 x0 = p0(:, 1);
 q0 = path_to_q(p0);
 r = norm(u);
 normal = x0*cos(vartheta) + ...
     (cross(x0, u)/norm(cross(x0, u)))*sin(vartheta);
 phi = r/sin(vartheta);
 x = x0*cos(phi) + cross(normal, x0)*sin(phi) +...
     normal*dot(normal, x0)*(1 - cos(phi));
 % the angle 'theta' that determines 'beta'
 theta = asin(sqrt(2/(1 + dot(x0, x)))*cos(vartheta));
 q = parallel_transport_PL(q0 + w, x0, x, theta);
 p = q_to_path(q, x);
end





