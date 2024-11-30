function power_turb = get_power(data,dist_from_turb,rho)
    power_turb = zeros(size(data.X_T,2),1);
    A = pi*data.R_T.^2;   % [m^2]

    for i = 1:size(data.X_T,2)
    % Interpolation to get the velocity plan at a dist from the turbine
    plane_in_front = interp3(data.z_u,data.y_u,data.x_u,data.u,data.z_u,data.y_u,data.X_T(i)-dist_from_turb);

    % Creation of a plan with Euclidian radius
    R_eucl = sqrt(data.z_u.^2 + data.y_u'.^2);

    % Keep only the value below the radius to obtain the circle aera before
    % the turbine
    mask = R_eucl <= data.R_T(i);
    plane_in_front(~mask) = NaN;

    % Power computation
    power_turb(i) = (0.5*rho*A(i)*nanmean(nanmean(plane_in_front))^3)./1e6;    % [MW]
    end
end