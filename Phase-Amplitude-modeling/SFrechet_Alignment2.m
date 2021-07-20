function [p, x, theta, paths, gammas, inv_gammas, SFF] = SFrechet_Alignment2(paths0, ...
    step1, step2, epsilon0, epsilon1, epsilon2, NN1, NN2, NN3)
%%% this algorithm finds the optimal 'p' that 
%%% is determined by the starting point 'x', which then 
%%% with the help of  'theta_Geodesics_Alignment' finds
%%% the geodesic baselines connecting 'x' with 'xj'
%%% which are determined by 'theta', where 
%%% with the help of 'Alignment' 
%%% optimally reparameterized paths 'paths', 
%%% respetively the reparameterizations 'gammas', and
%%% the iverse reparameterizations 'inv_gammas' 
%%% of original paths 'paths0' are derived
%%% thereby altogether defines TSRVC q at x
%%% so we get p = (x, q)
%%% INPUTS:
%%%% 'paths0' is the set of curves
%%%% 'epsilon0' is ideal stopping criteria for Alignment,
%%%% 'epsilon1' is ideal stopping criteria for theta_Geodesics_Alignment,
%%%% 'epsilon2' is ideal stopping criteria for gradient here,
%%%% 'step1' is moving step in theta_Geodesics_Alignment,
%%%% 'step2' is moving step here;
%%%%  the maximum number of iterations 'NN1',
%%%%  the maximum number of iterations when we don't perform
%%%%  steps 'NN2' and 'NN3';
%%% OUTPUTS;
%%%% 'p' is the sample Frechet mean,
%%%% 'x' is the starting point of 'p',
%%%% 'theta' is vector of parameters 'thetaj' 
%%%% that determines the baselines 'betaj' defining geodesics on 
%%%% the space of unparameterized curves,
%%%% 'paths' is set of reparameterized paths that are aligned to 'p',
%%%% 'gammas' is set of reparameterizations of 'paths0' to 'paths'
%%%% 'inv_gammas' is set of reparameterizations of 'paths' to 'paths0',
%%%% 'SFF' is the value of the sample Frechet Function. 
 N = length(paths0);% sample size 
 % initiate the starting point
 x = 0;
 for j=1:N
  pj = paths0{j};
  Xj = pj(:, 1);
  x = x + Xj;
 end
 if (x == 0)
    x = pj(:, 1);
 else
    x = x/norm(x);
 end
 % define the initial p-optimal baselines
 theta0 = zeros(1, N);
 step11 = step2;
 [theta, paths, gammas, inv_gammas] = theta_Geodesics_Alignment2(x, theta0,...
     paths0, step1, epsilon0, epsilon1, NN1/2, NN2/2, NN3/2);
 % find initial 'gradient' and the cost function 'SC'
 [gradient_x, SFF] = DxSampleFrechetFunction(x, theta, paths);
 disp(norm(gradient_x));
 % initialize counting variables
 N1 = 0;
 N2 = 0;
 N3 = 0;
 while (norm(gradient_x) > epsilon2)
  if (N1 < NN1)
   N1 = N1 + 1;   
   if (norm(step11*gradient_x) < pi)
    N3 = 0;   
    x1 = Exp_Sphere(x, - step11*gradient_x);
    SFF1 = SampleFrechetFunction(x1, theta, paths);
    if (SFF1 < SFF)
     N2 = 0;
     x = x1; 
     [theta, paths, gammas, inv_gammas] = theta_Geodesics_Alignment2(x, theta, ...
         paths0, step1, epsilon0, epsilon1, NN1/2, NN2/2, NN3/2);  
     [gradient_x, SFF] = DxSampleFrechetFunction(x, theta, paths);
     nn = norm(gradient_x);  
     X = sprintf('gradient with respect to x %e \n %d \n', nn, N1);
     disp(X);    
    else  
     if(N2 < NN2)   
     N2 = N2 + 1; step11 = step11/2;
     else
      nx = norm(gradient_x);
      X = sprintf('we performed more than N2=%d steps with norm of gradient with respect to x %e \n', NN2, nx);
      disp(X);
      break;
     end  
    end  
   else 
    if(N3 < NN3)
     N3 = N3 + 1;  step11 = step11/2;
    else
     nx = norm(gradient_x);
     X = sprintf('we performed more than N3=%d steps with norm of gradient with respect to x %e \n gradient is too big \n', NN3, nx);
     disp(X);
     break;
    end  
   end  
  else
   nx = norm(gradient_x);
   X = sprintf('we performed more than N1=%d steps with norm of gradient with respect to x %e \n', NN1, nx);
   disp(X);
   break;
  end 
 end
 p1 = paths0{1};
 q = 0*p1;
 for j=1:N   
  pj = paths{j};
  qj = path_to_q(pj);
  Xj = pj(:, 1);
  parqj = parallel_transport_PL(qj, Xj, x, theta(j));
  q = q + parqj;
 end
 q = q/N;
 p = q_to_path(q, x);
 SFF = 0;
 for j = 1:N
  p0j = paths0{j};   
  [thetaj, p0jn, pn, gamma, inv_gamma, square_length, gradient_theta] = theta_geodesic_alignment2(p0j, p, ...
       theta(j), step1, epsilon1, NN1/2, NN2/2, NN3/2);
  SFF = SFF + square_length;
  gammas{j} = gamma;
  inv_gammas{j} = inv_gamma;
  theta(j) = thetaj;
  paths{j} = p0jn;
 end    
 SFF = SFF/N;
end


