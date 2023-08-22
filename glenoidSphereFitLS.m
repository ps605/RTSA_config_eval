function [x_opt] = glenoidSphereFitLS(glenoid_stl, x0, glenoid_normal, glenoid_barycentre)

% Least squares iterative objective function
    function f=Objfun(x0)
        f=0; xe=x0(1); ye=x0(2); ze=x0(3);
        r=abs(x0(4));
        for i=1:n
            f=f+ 1*( (x(i)-xe)^2 + (y(i)-ye)^2 + (z(i)-ze)^2 - r^2)^2;
        end
    end

% Constraint for centre on lateral side of glenoid
    function [c, ceq] = centreCon(x0, glenoid_normal, glenoid_c)

       radius_vec = (x0(1:3) - glenoid_c)/norm(x0(1:3) - glenoid_c);
        
        ceq = [];
        %Inequality constraint for glenoid plane normal (lateral pointing)
        %and sphere radius vector to also point laterally 
        c(1) = -dot(glenoid_normal, radius_vec);
        % Inequality constraint to lateralise the centre of the LS sphere
        % (unstable after ~0.68) sphere retur radius r = 0
        c(2) = -(norm(x0(1:3) - glenoid_c) - 0.5*abs(x0(4)));
    end


n=size(glenoid_stl.Points,1); %number of data points

x = glenoid_stl.Points(:,1);
y = glenoid_stl.Points(:,2);
z = glenoid_stl.Points(:,3);

plot3(x,y,z,'o')
%least square fit to obtain unknowns x0,y0,z0 and r
%let the estimated variables be re,xe,ye and ze

options=optimset('MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',1e-6);
A = [];
b = [];
Aeq = [];%[0, 0, 0,1];
beq = [];%abs(x0(4));
lb = [];
ub = [];
nonlcon = @(x0)centreCon(x0, glenoid_normal, glenoid_barycentre);
x_opt=fmincon(@Objfun, ...
    x0, ...
    A,b,Aeq,beq,lb,ub,[nonlcon], ...
    options);

end