function [paths, gammas, inv_gammas] = Alignment(x, theta, paths0, epsilon0, NN1)
%%% this algorithm finds for the given starting point 'x' and
%%% the geodesic baselines connecting 'x' with 'xj',
%%% which are determined by 'theta',
%%% optimally reparameterized paths 'paths', 
%%% respetively the reparameterizations 'gammas', and
%%% the iverse reparameterizations 'inv_gammas' 
%%% of original paths 'paths0' are derived
%%% thereby altogether defines TSRVC q at x
%%% so we get p = (x, q)
%%% INPUTS:
%%%% 'x' is the starting point,
%%%% 'theta' is vector of parameters 'thetaj' that determine 
%%%% the geodesic baselines connecting 'x' with 'xj',
%%%% 'paths0' is the set of curves,
%%%% 'epsilon0' is ideal stopping criteria for here,
%%%%  the maximum number of iterations 'NN1';
%%% OUTPUTS;
%%%% 'paths' is set of reparameterized paths that are synchronized,
%%%% 'gammas' is set of reparameterizations of 'paths0' to 'paths'.
%%%% 'inv_gammas' is set of reparameterizations of 'paths' to 'paths0'. 
 N = length(paths0);% sample size
 SFF = SampleFrechetFunction(x, theta, paths0);
 paths = paths0;
 pj = paths0{1};
 q = 0*pj;
 for j=1:N
  pj = paths0{j};
  Xj = pj(:, 1);
  qj = path_to_q(pj);
  par_qj = parallel_transport_PL(qj, Xj, x, theta(j));
  q = q + par_qj;
 end
 q0 = q/N;
 T = size(pj, 2);
 t = linspace(0, 1, T);
 q = 0*q0;
 for j=1:N
  p0j = paths0{j};
  Xj = p0j(:, 1);
  q0j = path_to_q(p0j);
  par_q0j = parallel_transport_PL(q0j, Xj, x, theta(j));
  [G, T1] = DynamicProgrammingQ2(q0, t, par_q0j, t, t, t, 0);
  gamma = interp1(T1, G, t);
  p0jn = Group_Action_by_Gamma_p(p0j, gamma);
  paths1{j} = p0jn;
  gammas1{j} = gamma;
  q0jn = path_to_q(p0jn);
  par_q0jn = parallel_transport_PL(q0jn, Xj, x, theta(j));
  q = q + par_q0jn;
  inv_gamma = invertGamma(gamma);
  inv_gamma = (inv_gamma-inv_gamma(1))/(inv_gamma(end)-inv_gamma(1));
  inv_gammas1{j} = inv_gamma;
 end 
 SFF1 = SampleFrechetFunction(x, theta, paths1);
 maxgamma = SFF - SFF1;
 if (maxgamma < epsilon0)
  for j=1:N
   pj = paths0{j};
   paths{j} = pj;
   gammas{j} = t;
   inv_gammas{j} = t;
  end    
 else 
  SFF = SFF1;   
  q0 = q/N;
  N1 = 0;
  paths = paths1;
  gammas = gammas1;
  inv_gammas = inv_gammas1;
  while(maxgamma > epsilon0) 
   if (N1 < NN1)
    N1 = N1+1;   
    q = 0*q0;
    for j=1:N
     p0j = paths0{j};
     Xj = pj(:, 1);
     q0j = path_to_q(p0j);
     par_q0j = parallel_transport_PL(q0j, Xj, x, theta(j));
     [G, T1] = DynamicProgrammingQ2(q0, t, par_q0j, t, t, t, 0);
     gamma = interp1(T1, G, t);
     p0jn = Group_Action_by_Gamma_p(p0j, gamma);
     paths1{j} = p0jn;
     gammas1{j} = gamma;
     q0jn = path_to_q(p0jn);
     par_q0jn = parallel_transport_PL(q0jn, Xj, x, theta(j));
     q = q + par_q0jn;
     inv_gamma = invertGamma(gamma);
     inv_gamma = (inv_gamma-inv_gamma(1))/(inv_gamma(end)-inv_gamma(1));
     inv_gammas1{j} = inv_gamma;
    end
    SFF1 = SampleFrechetFunction(x, theta, paths1);
    if (SFF1 < SFF)
     maxgamma = SFF - SFF1;
     q0 = q/N;
     paths = paths1;
     gammas = gammas1;
     inv_gammas = inv_gammas1;
     SFF = SFF1;
     nn = maxgamma;  
     X = sprintf('difference of Cost functions %e \n at iteration %d \n', nn, N1);
     disp(X);
    else
     nn = maxgamma;   
     X = sprintf('Cost function is not reduced \n at iteration %d \n with last difference %e \n', N1, nn);
     disp(X);
     break;
    end 
   else
    nx = maxgamma;
    X = sprintf('we performed more than %d step \n with difference of Cost functions %e \n', NN1, nx);    
    disp(X);
    break;
   end 
  end
 end
end