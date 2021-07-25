function [theta, p1n, p2n, gamma, inv_gamma, square_length, gradient_theta] = theta_geodesic_alignment2(p1, p2, ...
   theta0, step, epsilon, NN1, NN2, NN3)
%%% the function is finding optimal 'theta' 
%%% that determines the the p-optimal baseline which generates
%%% the geodesic between p1 and p2 and the space of unparameterized curves
%%% the inputs: two curves 'p1' and 'p2', the initial parameter 'theta0',
%%%             the degree of accurace 'epsilon', the step size 'step',
%%%             the maximum number of iterations 'NN1',
%%%             the maximum number of iterations when we don't perform
%%%             steps 'NN2' and 'NN3'
%%% the outputs: the 'theta' that determines baseline 'beta',
%%%              the optimal square length 'square_length',
%%%              the gradient of square length 'gradient_theta'
%%%              with respect to 'theta'
%%%              the reparameterization 'gamman' for 'p1', 
%%%              the optimally reparameterized curve 'p1n',
%%%              the reparameterization 'inv_gamman' for 'p2', 
%%%              the optimally reparameterized curve 'p2n',
% we have starting points on S2
 T = size(p1, 2);
 t = linspace(0, 1, T);
 x1 = p1(:, 1);
 x2 = p2(:, 1);
 if ( norm(x1 - x2) < 0.00000001 )
  theta = 0;
  p1n = p1; p2n = p2; 
  gamma = t; inv_gamma = t;
  square_length = SquareLength(p1, p2, theta); gradient_theta = 0;
  [square_length_next, p1n_next, p2n_next, gamma_next, inv_gamma_next] = SquareLength_Alignment(p1, p2, theta);
  if(square_length_next < square_length)
   square_length = square_length_next;
   p1n = p1n_next;
   p2n = p2n_next;
   gamma = gamma_next;
   inv_gamma = inv_gamma_next;  
  end
 else 
  theta = theta0;
  p1n = p1; p2n = p2; 
  gamma = t; inv_gamma = t;
  [gradient_theta, square_length] = DthetaSquareLength(p1, p2, theta);  
  [gradient_theta_next, p1n_next, p2n_next, gamma_next, inv_gamma_next, square_length_next] = DthetaSquareLength_Alignment(p1, ...
      p2, theta0);
  if(square_length_next < square_length)
   square_length = square_length_next;
   p1n = p1n_next;
   p2n = p2n_next;
   gamma = gamma_next;
   inv_gamma = inv_gamma_next;  
   gradient_theta = gradient_theta_next;
  end
  if (abs(gradient_theta) <= epsilon)
   theta = theta0;
   %llength = SquareLength(p1n, p2, theta);
  else
   theta = theta0;
   M = 8;
   t = linspace(-pi/2, pi/2, M);
   for k=1:M
    sl(k) =  SquareLength(p1, p2, t(k));
    [square_length_next2(k), p1n_next2, p2n_next2, gamma_next2, inv_gamma_next2] = SquareLength_Alignment(p1, p2, t(k));
   end
   [sl2 , mink] = min(square_length_next2);
   theta2 = t(mink);
   if (square_length > sl2)
    theta = theta2;
    [gradient_theta, p1n, p2n, gamma, inv_gamma, square_length] = DthetaSquareLength_Alignment(p1, ...
      p2, theta);
   end 
   % set counting parameters
   N1 = 0;
   while (abs(gradient_theta)>epsilon) 
    if (N1 <= NN1)
     N1 = N1 + 1;  
     if (square_length == SquareLength(p1n, p2, theta))
      [theta, square_length_next, gradient_theta] = theta_geodesic(p1n, p2, theta, step, epsilon, NN1*3/5, NN2*3/5, NN3*3/5);
     else
      [theta, square_length_next, gradient_theta] = theta_geodesic(p1, p2n, theta, step, epsilon, NN1*3/5, NN2*3/5, NN3*3/5);  
     end
     if (square_length_next < square_length)
      square_length = square_length_next;
     else
       break;
     end  
     [gradient_theta_next, p1n_next, p2n_next, gamma_next, inv_gamma_next, square_length_next] = DthetaSquareLength_Alignment(p1, ...
         p2, theta);
     if (square_length > square_length_next)
      square_length = square_length_next;
      p1n = p1n_next;
      p2n = p2n_next;
      gamma = gamma_next;
      inv_gamma = inv_gamma_next;
      gradient_theta = gradient_theta_next;
      X = sprintf('the gradient with respect to theta \n with alignment  %e \n', abs(gradient_theta));
      disp(X);
     end
    else
     X = sprintf('we reached the number of iterations %d with alignment \n and with gradient with respect to theta %e \n', NN1, abs(gradient_theta));
     disp(X);
     break
    end
   end
  end
 end 
end
