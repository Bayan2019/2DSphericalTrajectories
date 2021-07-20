function [gradient_theta, p1n, p2n, gamma, inv_gamma, square_length] = DthetaSquareLength_Alignment(p1, p2, theta)
%%% the derivative 'gradient_theta' of the square length 'SquareLength_Alignment'
%%% with respect to 'theta' on the space of unparameterized curves
%%% Inputs: two curves 'p1' and 'p2', 'theta' which determines 
%%%         the baselines 'beta'
%%% Outputs: the derivative 'gradient_theta' 
%%%          of the square length 'SquareLength_Alignment'
%%%          with respect to 'theta'
%%%          and the square length 'square_length'
%%%          the reparameterization 'gamman' for 'p1', 
%%%          the optimally reparameterized curve 'p1n',
%%%          the reparameterization 'inv_gamman' for 'p2', 
%%%          the optimally reparameterized curve 'p2n',
% the starting points
 x1 = p1(:, 1);
 x2 = p2(:, 1);
 % and SRVC at tangent spaces
 q1 = path_to_q(p1);
 q2 = path_to_q(p2);
 % the time grid
 T = size(q1, 2);
 t = linspace(0, 1, T);
 if ( norm(x1 - x2) < 0.000001 )
  gradient_theta = 0;
  % alignment of qi
  [G, T1] = DynamicProgrammingQ2(q2, t, q1, t, t, t, 0);
  gamma = interp1(T1, G, t);
  p1n = Group_Action_by_Gamma_p(p1, gamma);
  inv_gamma = invertGamma(gamma);
  inv_gamma = (inv_gamma-inv_gamma(1))/(inv_gamma(end)-inv_gamma(1));
  p2n = Group_Action_by_Gamma_p(p2, inv_gamma);
  % computing the length
  square_length1n = SquareLength(p1n, p2, theta);
  square_length2n = SquareLength(p1, p2n, theta);
  square_length = min(square_length1n, square_length2n);
 else
  % parallel transport
  par_q1 = parallel_transport_PL(q1, x1, x2, theta);
  % alignment of the parallel transport parqi
  [G, T1] = DynamicProgrammingQ2(q2, t, par_q1, t, t, t, 0);
  gamma = interp1(T1, G, t);
  p1n = Group_Action_by_Gamma_p(p1, gamma);
  inv_gamma = invertGamma(gamma);
  p2n = Group_Action_by_Gamma_p(p2, inv_gamma);
  % gradient
  square_length1n = SquareLength(p1n, p2, theta);
  square_length2n = SquareLength(p1, p2n, theta);
  if (square_length1n <= square_length2n)
   [gradient_theta, square_length] = DthetaSquareLength(p1n, p2, theta);
  else
   [gradient_theta, square_length] = DthetaSquareLength(p1, p2n, theta);  
  end    
 end
end