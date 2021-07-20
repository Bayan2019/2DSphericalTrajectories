function [theta, SFF] = theta_Geodesics(x, theta0, paths, step, epsilon, NN1, NN2, NN3)
%%% this algorithm finds for the given starting point 'x'
%%% the geodesic baselines connecting 'x' with 'xj'
%%% which are determined by 'theta', where 
%%% thereby altogether defines TSRVC q at x
%%% so we get p = (x, q)
%%% INPUTS:
%%%% 'x' is the starting point,
%%%% 'theta0' is the initial values, 
%%%% 'paths' is the set of curves,
%%%% 'epsilon1' is ideal stopping criteria here,
%%%% 'step1' is moving step in here,
%%%%  the maximum number of iterations 'NN1',
%%%%  the maximum number of iterations when we don't perform
%%%%  steps 'NN2' and 'NN3'
%%% OUTPUTS;
%%%% 'theta' is vector of parameters 'thetaj' 
%%%%         that determine the baselines 'betaj' defining geodesics on 
%%%%         the space of parameterized curves,
%%%% 'SFF' is the value of the sample Frechet Function. 
 [gradient_theta, SFF] = DthetaSampleFrechetFunction(x, theta0, paths);
 step1 = step;
 theta = theta0;
 N1 = 0;
 N2 = 0;
 N3 = 0; nn = 1;
 while (norm(gradient_theta)>epsilon)
  if (N1 < NN1)
   N1 = N1 + 1;
   theta1 = theta - step1*gradient_theta;
   if (all(theta1 >= -pi/2) && all(theta1 <= pi/2))
    N3 = 0;
    SFF1 = SampleFrechetFunction(x, theta1, paths);
    if (SFF1 < SFF)
     N2 = 0;
     theta = theta1;  
     [gradient_theta, SFF] = DthetaSampleFrechetFunction(x, theta, paths);
     %nn = norm(step1*gradient_theta);  
     %X = sprintf('gradient with respect to theta without alignment %e \n at the iteration %d \n', ...
     %            norm(gradient_theta), N1);
     %disp(X);
    else
     if(N2 < NN2)
      N2 = N2 + 1; nn = norm(step1*gradient_theta);
      step1 = step1/2;  
     else
      nn = norm(step1*gradient_theta);   
      X = sprintf('we performed more than %d steps without moving \n and get gradient with respect to theta without alignment %e \n at the neighborhood %e \n', ...
          NN2, norm(gradient_theta), nn);
      disp(X);
      break   
     end    
    end    
   else
    if(N3 < NN3)
     N3 = N3 + 1;
     step1 = step1/2;
    else  
     X = sprintf('we performed more than %d steps and \n get gradient with respect to theta without alignment %e around -pi/2 or pi/2 \n', ...
         NN3, norm(gradient_theta));
     disp(X);
     break  
    end    
   end
  else
   X = sprintf('we performed more than %d steps and gradient with respect to theta without alignment %e \n at the neighborhood %e \n', ...
       NN1, norm(gradient_theta), nn);
   disp(X);
   break
  end
 end
end
