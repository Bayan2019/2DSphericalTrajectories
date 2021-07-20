function [theta, paths, gammas, inv_gammas] = theta_Geodesics_Alignment2(x, theta0, paths0, step, epsilon0, epsilon1, NN1, NN2, NN3)
%%% this algorithm finds for the given starting point 'x'
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
%%%% 'x' is the starting point,
%%%% 'theta0' is the initial values, 
%%%% 'paths0' is the set of curves,
%%%% 'epsilon0' is ideal stopping criteria for Alignment,
%%%% 'epsilon1' is ideal stopping criteria here,
%%%% 'step1' is moving step in here,
%%%%  the maximum number of iterations 'NN1',
%%%%  the maximum number of iterations when we don't perform
%%%%  steps 'NN2' and 'NN3'
%%% OUTPUTS;
%%%% 'theta' is vector of parameters 'thetaj' 
%%%% that determine the baselines 'betaj' defining geodesics on 
%%%% the space of unparameterized curves,
%%%% 'paths' is set of reparameterized paths that are synchronized,
%%%% 'gammas' is set of reparameterizations of 'paths0' to 'paths'.
%%%% 'inv_gammas' is set of reparameterizations of 'paths' to 'paths0'. 
 [paths, gammas, inv_gammas] = Alignment(x, theta0, paths0, epsilon0, NN1);
 [gradient_theta, SFF] = DthetaSampleFrechetFunction(x, theta0, paths);
 step1 = step;
 theta = theta0;
 N1 = 0;
 N2 = 0;
 while (norm(gradient_theta) > epsilon1)
  if (N1 < NN1)
   N1 = N1+1;
   [theta1, SFF1] = theta_Geodesics(x, theta, paths, step, epsilon1, NN1, NN2, NN3);
   if (SFF1 < SFF)
    N2 = 0;
    theta = theta1;
    [paths, gammas, inv_gammas] = Alignment(x, theta, paths0, epsilon0, NN1);        
    [gradient_theta, SFF] = DthetaSampleFrechetFunction(x, theta, paths);
    nn = norm(gradient_theta);  
    X = sprintf('gradient with respect to theta %e \n at the iteration %d \n', nn, N1);
    disp(X);   
   else      
    if(N2 < NN2)
     N2 = N2 + 1; step1 = step1/2;
    else
     nn = norm(gradient_theta);
     X = sprintf('we performed more than N2=%d steps and get gradient with respect to theta %e \n', NN2, nn);
     disp(X);
     break
    end          
   end  
  else
   nn = norm(gradient_theta);
   X = sprintf('we performed more than N1=%d steps and gradient with respect to theta %e \n', NN1, nn);
   disp(X);
   break
  end 
 end
end
