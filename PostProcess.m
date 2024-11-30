clear all; close all; clc;
%% Post-processing file 
rho = 1.225;    % Air density [kg/m^3]

%% Load of the data from Python script
data4 = load('simulation_1466_4ms.mat');
data5 = load('simulation_1469_5ms.mat');
data6 = load('simulation_1471_6ms.mat');
data7 = load('simulation_1473_7ms.mat');
data = load('simulation_1475_8ms.mat');     % Ref case for the post processing

% Convert each field to double
fieldNames = fieldnames(data4);
for i = 1:numel(fieldNames)
    data4.(fieldNames{i}) = double(data4.(fieldNames{i}));
end

% Convert each field to double
fieldNames = fieldnames(data5);
for i = 1:numel(fieldNames)
    data5.(fieldNames{i}) = double(data5.(fieldNames{i}));
end

% Convert each field to double
fieldNames = fieldnames(data6);
for i = 1:numel(fieldNames)
    data6.(fieldNames{i}) = double(data6.(fieldNames{i}));
end

% Convert each field to double
fieldNames = fieldnames(data7);
for i = 1:numel(fieldNames)
    data7.(fieldNames{i}) = double(data7.(fieldNames{i}));
end

% Convert each field to double
fieldNames = fieldnames(data);
for i = 1:numel(fieldNames)
    data.(fieldNames{i}) = double(data.(fieldNames{i}));
end


x = data.x_p; y = data.y_p; z = data.z_p;

%% Visual plot of 2D mesh and speed countour
% Create a 3D grid using meshgrid
[X, Z] = meshgrid(x, z);

% Plot the mesh in 2D space
figure(1)
subplot(211)
hold on;

% Draw grid lines parallel to the x-axis
for i = 1:length(z)
    plot(x, z(i) * ones(size(x)), 'k');  % Horizontal lines
end

% Draw grid lines parallel to the y-axis
for i = 1:length(x)
    plot(x(i) * ones(size(z)), z, 'k');  % Vertical lines
end

% Turbine representation
for i = 1:size(data.X_T,2)
    turb = linspace(-data.R_T(i),data.R_T(i),100);
    plot(data.X_T(i)*ones(size(turb)),turb,'-r','LineWidth',2)
end

xlabel('X [m]'); ylabel('Z [m]');
title('2D Mesh [y = 0 plan]');
grid on; box on;


% Plot of the contour
subplot(212)
hold on; 

contourf(data.x_u,data.z_u,squeeze(data.u(:,72,:)))

% Turbine representation
for i = 1:size(data.X_T,2)
    turb = linspace(-data.R_T(i),data.R_T(i),100);
    plot(data.X_T(i)*ones(size(turb)),turb,'-r','LineWidth',2)
end

xlabel('X [m]'); ylabel('Z [m]');
title('Velocity Countour [y = 0 plan] [m/s]');
grid on; box on;
hColorbar = colorbar('east');


%% Plot of centerline u and p
u_centerline = zeros(size(data.x_u,2),size(data.X_T,2));
p_centerline = zeros(size(data.x_p,2),size(data.X_T,2));

% Interpolation over Y and Z to have the hub height value
for i = 1:size(data.X_T,2)
    u_centerline(:,i) = interp3(data.z_u,data.y_u,data.x_u,data.u,data.Z_T(i),data.Y_T(i),data.x_u);
    p_centerline(:,i) = interp3(data.z_p,data.y_p,data.x_p,data.p,data.Z_T(i),data.Y_T(i),data.x_p);
end

% Plot of the centerline value for each turbine (useful if there is a
% different Y or Z than o)
figure(2)
subplot(211)
plot(data.x_u,u_centerline(:,1))   
xlabel('X [m]'); ylabel('HH wind speed [m/s]');
title('Centerline streamwise velocity');
grid on; box on; hold on; axis tight;

subplot(212)
plot(data.x_p,p_centerline(:,1)) 
xlabel('X [m]'); ylabel('Pressure [?]');
title('Centerline Pressure');
grid on; box on; hold on; axis tight;

%% Power Conversion Table (to get rid of induction)
dist_from_turb = 7*data.R_T;    % [m]
power_at_rotor = [get_power(data4,dist_from_turb,rho) get_power(data5,dist_from_turb,rho) get_power(data6,dist_from_turb,rho) get_power(data7,dist_from_turb,rho) get_power(data,dist_from_turb,rho)];
power_bf_rotor = [get_power(data4,0,rho) get_power(data5,0,rho) get_power(data6,0,rho) get_power(data7,0,rho) get_power(data,0,rho)];

figure(3)
plot(power_bf_rotor,power_at_rotor,'or','LineWidth',2)

ylabel(sprintf('Power at %2.1f D [MW]',dist_from_turb/(2*data.R_T))); xlabel('Power at rotor [MW]');
title('Power regression fit');
grid on; box on; hold on; axis tight;

% Get the linear reg fit
p = polyfit(power_bf_rotor,power_at_rotor,1);

% plot of linear fit
yfit = p(1).*linspace(power_bf_rotor(1),power_bf_rotor(end),100) + p(2);
plot(linspace(power_bf_rotor(1),power_bf_rotor(end),100),yfit,'b--','LineWidth',1.5)



%% Power computation (Area average)
power_turb = get_power(data,dist_from_turb,rho);

% Computation of normalized power
power_norm = 1e6*power_turb./(0.5.*rho.*(pi*data.R_T.^2)'.*data.U_inf.^3);




