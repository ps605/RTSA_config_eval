function [c, ceq] = sphere_func_con(p_sim, p_sim_0, radius)

c = (p_sim(1) - p_sim_0(1))^2 + (p_sim(2) - p_sim_0(2))^2 + (p_sim(3) - p_sim_0(3))^2 - radius^2;

ceq =[]; 

% % figure(103);
% % title('c')
% % scatter(1, c, 'o', 'filled', 'cyan')
% % % ylim([-0.05 0.05])
% % hold on
end