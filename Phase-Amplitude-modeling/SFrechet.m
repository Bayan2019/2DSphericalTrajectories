function [p, x, theta, SFF] = SFrechet(paths0, step1, step2, epsilon1, epsilon2, NN1, NN2, NN3)
%%% this algorithm finds the optimal 'p' that 
%%% is determined by the starting point 'x', which then 
%%% with the help of  'theta_Geodesics' finds
%%% the geodesic baselines connecting 'x' with 'xj'
%%% which are determined by 'theta',
%%% thereby altogether defines TSRVC q at x
%%% so we get p = (x, q)
%%% INPUTS:
%%%% 'paths0' is the set of curves
%%%% 'epsilon1' is ideal stopping criteria for theta_Geodesics,
%%%% 'epsilon2' is ideal stopping criteria for gradient here,
%%%% 'step1' is moving step in theta_Geodesics,
%%%% 'step2' is moving step here;
%%%%  the maximum number of iterations 'NN1',
%%%%  the maximum number of iterations when we don't perform
%%%%  steps 'NN2' and 'NN3';
%%% OUTPUTS;
%%%% 'p' is the sample Frechet mean,
%%%% 'x' is the starting point of 'p',
%%%% 'theta' is vector of parameters 'thetaj' 
%%%%         that determines the baselines 'betaj' defining geodesics on 
%%%%         the space of parameterized curves,
%%%% 'SFF' is the value of the Sample Frechet Function. 
 N = length(paths0); % sample size
 x0 = 0;
 for j=1:N
  p0j = paths0{j};
  Xj = p0j(:, 1);
  x0 = x0 + Xj;
 end
 if (x0 == 0)
  x0 = p0j(:, 1);
 else
  x0 = x0/norm(x0);
 end
 theta0 = zeros(1, N);
 step11 = step2;
 [theta, SFF] = theta_Geodesics(x0, theta0, paths0, step1, epsilon1,  NN1/2, NN2/2, NN3/2);
 [gradient_x, SFF] = DxSampleFrechetFunction(x0, theta, paths0);
 x = x0;
 disp(norm(gradient_x));
 % initiate counting variables
 N1 = 0;
 N2 = 0;
 N3 = 0;
 while (norm(gradient_x) > epsilon2)
  if (N1 < NN1)
   N1 = N1 + 1;   
   if (norm(step11*gradient_x) < pi)
    N3 = 0;   
    x1 = Exp_Sphere(x, - step11*gradient_x);
    SFF1 = SampleFrechetFunction(x1, theta, paths0);
    if (SFF1 < SFF)
     N2 = 0; 
     x = x1;  
     [theta, SFF] = theta_Geodesics(x, theta, paths0, step1, epsilon1, NN1/2, NN2/2, NN3/2);
     [gradient_x, SFF] = DxSampleFrechetFunction(x, theta, paths0);
     nx = norm(gradient_x);   
     X = sprintf('gradient with respect to x without alignment %e \n', nx);
     disp(X);
    else  
     if(N2 < NN2)   
      N2 = N2 + 1; step11 = step11/2;
     else
      nx = norm(gradient_x);
      X = sprintf('we performed more than N2=%d step with norm of gradient with respect to x without alignment %e \n', NN2, nx);
      disp(X);
      break
     end
    end
   else
    if(N3 < NN3)
     N3 = N3 + 1; step11 = step11/2;
    else
     nx = norm(gradient_x);
     X = sprintf('we performed more than N3=%d step with norm of gradient with respect to x without alignment %e \n around -pi/1 and pi/2 \n', NN3, nx);
     disp(X);
     break;
    end
   end
  else
   nx = norm(gradient_x);
   X = sprintf('we performed more than N1=%d step with norm of gradient with respect to x without alignment %e \n', NN1, nx);
   disp(X);
   break;
  end
 end
 q = 0*p0j;
 for j=1:N   
  p0j = paths0{j};
  q0j = path_to_q(p0j);
  Xj = p0j(:, 1);
  parq0j = parallel_transport_PL(q0j, Xj, x, theta(j));
  q = q + parq0j;
 end
 q = q/N;
 p = q_to_path(q, x);
 SFF = 0;
 for j = 1:N
  p0j = paths0{j};   
  [thetaj, square_length, gradient_theta] = theta_geodesic(p0j, p, ...
       theta(j), step1, epsilon1, NN1/2, NN2/2, NN3/2);
  SFF = SFF + square_length;
  theta(j) = thetaj;
 end 
 SFF = SFF/N;
end
