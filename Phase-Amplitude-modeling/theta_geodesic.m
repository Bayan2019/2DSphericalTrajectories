function [theta, square_length, gradient_theta] = theta_geodesic(p1, ...
    p2, theta0, step, epsilon, NN1, NN2, NN3)
%%% the function is finding 'theta' that determines the baseline which generates
%%% the geodesic between p1 and p2 and the space of parameterized curves
%%% the inputs: two curves 'p1' and 'p2', the initial parameter '\theta',
%%%             the degree of accurace 'epsilon', the step size 'step',
%%%             the maximum number of iterations 'NN1',
%%%             the maximum number of iterations when we don't perform
%%%             steps 'NN2' and 'NN3'
%%% the outputs: the 'theta' that determines baseline 'beta',
%%%              the optimal square length 'square_length',
%%%              the gradient of square length 'gradient_theta'
%%%              with respect to 'theta' 
% we have starting points on S2
x1 = p1(:, 1);
x2 = p2(:, 1);
% if we assume that q1 and q2 have the same length 
% or determined at the same finite set of points in [0, 1]
 if ( norm(x1 - x2) < 0.0000001 )
  theta = 0;
  square_length = SquareLength(p1, p2, theta);
  gradient_theta = 0;
 else    
  [gradient_theta, square_length] = DthetaSquareLength(p1, p2, theta0);
  if (abs(gradient_theta) <= epsilon)
   theta = theta0;
  else
   step1 = step;
   theta = theta0;
   N1 = 0;
   N2 = 0;
   N3 = 0; nn = 1;
   while (abs(gradient_theta) > epsilon)
    if (N1 <= NN1)    
     theta1 = theta - step1*gradient_theta; 
     N1 = N1 + 1;
     if ((theta1 >= -pi/2)&&(theta1 <= pi/2))
      N3 = 0;   
      square_length1 = SquareLength(p1, p2, theta1);       
      if (square_length1 < square_length)
       N2 = 0;    
       theta = theta1; 
       [gradient_theta, square_length] = DthetaSquareLength(p1, p2, theta1);
       % X = sprintf('the gradient with respect to theta %e \n at the iteration %d \n without an alignment\n', abs(gradient), N1);
       % disp(X);
      else
       if(N2 > NN2)
        nn = abs(step1*gradient_theta);
        X = sprintf('we reached the number of iterations N2 = %d, \n we are very close in neighborhood of %e \n with the gradient %e \n without alignment \n', ...
            NN2, nn, abs(gradient_theta));
        disp(X);
        break
       else
        N2 = N2 + 1; nn = abs(step1*gradient_theta);
        step1 = step1/2;
       end
      end
     else
      if (N3 > NN3)
       X = sprintf('we reached the number of iterations N3 = %d, \n we reached the bound -pi/2 or pi/2 \n with the gradient %e \n without alignment \n', ...
           NN3, abs(gradient_theta));
       disp(X);
       break
      else
       N3 = N3 + 1;
       step1 = step1/2;
      end
     end
    else 
     X = sprintf('we reached the number of iterations %d \n with gradient %e \n without alignment \n a the neighborhood %e \n', ...
         NN1, abs(gradient_theta), nn);
     disp(X);
     break
    end     
   end
  end
 end
end
