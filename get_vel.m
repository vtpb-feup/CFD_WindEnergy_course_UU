function vel_disk = get_vel(data,dist_from_turb,X_T,R_T)
    vel_disk = zeros(size(X_T,2),1);

    for i = 1:size(X_T,2)
    % Interpolation to get the velocity plan at a dist from the turbine
    plane_in_front = interp3(data.y_u,data.z_u,data.x_u,data.u,data.z_u,data.y_u,X_T(i)-dist_from_turb);

    % Creation of a plan with Euclidian radius
    R_eucl = sqrt(data.z_u.^2 + data.y_u'.^2);

    % Keep only the value below the radius to obtain the circle aera before
    % the turbine
    mask = R_eucl <= R_T(i);
    plane_in_front(~mask) = NaN;

    % Power computation
    vel_disk(i) = nanmean(nanmean(plane_in_front)); % [MW]
    end
end