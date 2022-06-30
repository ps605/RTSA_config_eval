function [c, ceq] =plane_func_d(y_m, plane_parameters)

c = [];
ceq = plane_parameters(1)*y_m(1) + plane_parameters(2)*y_m(2) + plane_parameters(3)*y_m(3) + plane_parameters(4);

