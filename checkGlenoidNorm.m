function [glenoid_normal, stl_scap] = checkGlenoidNorm(x, y, z, glenoid_normal, glenoid_barycentre)
%% Check which way the Plane norm (Z axis) is pointing
% Get the barycentre of glenoid verteces
x_sca_mean = mean(x,'all');
y_sca_mean = mean(y,'all');
z_sca_mean = mean(z,'all');

scap_barycentre = [x_sca_mean y_sca_mean z_sca_mean];

scatter3(scap_barycentre(1), scap_barycentre(2), scap_barycentre(3),'cyan', 'filled','o')

% Pass .stl x-y-z to structure to then export
stl_scap = struct('x', x,...
    'y', y,...
    'z', z);

% What is the angle between the norm and the vector from the norm origine to glenoid vertex barycentre

% Vector from glenoid barycentre to humeral verteces barrycentre
gle_bary_vec = scap_barycentre - glenoid_barycentre;

% Connect with line to visualise vector
line([glenoid_barycentre(1) scap_barycentre(1)],...
    [glenoid_barycentre(2) scap_barycentre(2)],...
    [glenoid_barycentre(3) scap_barycentre(3)], ...
    'LineWidth',4,'Color','cyan');

norm_check_angle = vrrotvec(gle_bary_vec, glenoid_normal);

% If angle between vectors is 0<gonia<=90 flip normal about plane
% This is the Z-axis
if norm_check_angle(4) >= 0 && norm_check_angle(4) < pi/2

    % Negate Z normal to flip about plane to point out of the glenoid
    glenoid_normal = - glenoid_normal;

end
end