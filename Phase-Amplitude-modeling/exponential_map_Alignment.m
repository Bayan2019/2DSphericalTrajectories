function [p, theta, vartheta, gradient_theta2] = exponential_map_Alignment(p0, u, w, step , epsilon, NN1, NN2, NN3)
 T = size(p0, 2);
 t = linspace(0, 1, T);
 q0 = path_to_q(p0);
 wq0 = w + q0;
 [G0, T0] = DynamicProgrammingQ2(q0, t , wq0, t, t, t, 0);
 gam = interp1(T0, G0, t);
 wq0_n = Group_Action_by_Gamma_q(wq0, gam);
 w_n = wq0_n - q0;
 [p, theta, vartheta, gradient_theta2] = exponential_map(p0, u, w_n, step , epsilon, NN1, NN2, NN3);
end

