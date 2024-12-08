clear all; close all; clc;
%% Post-processing file
% Adapt these values in function of the file your are post-processing
rho = 1.225;    % Air density [kg/m^3]
X_T = [1000 1300 1600 1900];
Y_T = [0 0 0 0];
Z_T = [0 0 0 0];
R_T = [50 50 50 50];
U_inf = 8;

%% Load of the data from Python script
data = load('simulation_R4_Sx3.mat');     % Ref case for the post processing

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
for i = 1:size(X_T,2)
    turb = linspace(-R_T(i),R_T(i),100);
    plot(X_T(i)*ones(size(turb)),turb,'-r','LineWidth',2)
end

xlabel('X [m]'); ylabel('Z [m]');
title('2D Mesh [y = 0 plan]');
grid on; box on;


% Plot of the contour
subplot(212)
hold on; 

contourf(data.x_u,data.z_u,squeeze(data.u(:,size(data.u,2)/2,:)))

% Turbine representation
for i = 1:size(X_T,2)
    turb = linspace(-R_T(i),R_T(i),100);
    plot(X_T(i)*ones(size(turb)),turb,'-r','LineWidth',2)
end

xlabel('X [m]'); ylabel('Z [m]');
title('Velocity Countour [y = 0 plan] [m/s]');
grid on; box on;
hColorbar = colorbar('east');


%% Plot of centerline u and p
u_centerline = zeros(size(data.x_u,2),size(X_T,2));
p_centerline = zeros(size(data.x_p,2),size(X_T,2));

% Interpolation over Y and Z to have the hub height value
for i = 1:size(X_T,2)
    u_centerline(:,i) = interp3(data.y_u,data.z_u,data.x_u,data.u,Z_T(i),Y_T(i),data.x_u);
    p_centerline(:,i) = interp3(data.y_p,data.z_p,data.x_p,data.p,Z_T(i),Y_T(i),data.x_p);
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


%% Power computation (Area average)
power_turb = (0.5.*rho.*(pi*R_T.^2)'.*get_vel(data,0,X_T,R_T).^3)/1e6;

%% Computation of Cp

% computation of Ref wind speed
p = [0.7024 0.1363];    % Linear coefficient from previous sim

U_ref = get_vel(data,0,X_T,R_T)/p(1) - p(2); 
Cp = 1e6*power_turb./(0.5.*rho.*(pi*R_T.^2)'.*U_ref.^3);
Cp(1) = 1e6*power_turb(1)./(0.5.*rho.*(pi*R_T(1).^2)'.*U_inf.^3);

%% Blockage analysis
% Import of wind speed of first turbine
WindSpeed_FirstTurbine = load("WindSpeed_firstTurbine.mat");

figure(3)
subplot(121)
Sx = 3*2*R_T(1); Sy = 4*2*R_T(1);
plot([1 2 4],[WindSpeed_FirstTurbine.U_T1_inf WindSpeed_FirstTurbine.U_R2_Sx3 WindSpeed_FirstTurbine.U_R4_Sx3]./WindSpeed_FirstTurbine.U_T1_inf,'ro--')
xlabel('N_{rows}'); ylabel('U/U_{ref}'); ylim([0.9905 1])
grid on; box on; hold on;
title('S_{x} = 3D')
Segalini_law = (1-0.097*((Sx*Sy/(2*R_T(1))^2)^-0.9)*(1-exp(0.88-0.88.*linspace(1,4,100))));
plot(linspace(1,4,100),Segalini_law,'b')


subplot(122)
Sx = 6*2*R_T(1); Sy = 4*2*R_T(1); 
plot([1 2 4],[WindSpeed_FirstTurbine.U_T1_inf WindSpeed_FirstTurbine.U_R2_Sx6 WindSpeed_FirstTurbine.U_R4_Sx6]./WindSpeed_FirstTurbine.U_T1_inf,'ro--')
xlabel('N_{rows}'); ylim([0.9905 1]);
grid on; box on; hold on;
title('S_{x} = 6D')
Segalini_law = (1-0.097*((Sx*Sy/(2*R_T(1))^2)^-0.9)*(1-exp(0.88-0.88.*linspace(1,4,100))));
plot(linspace(1,4,100),Segalini_law,'b')
legend('CFD Data','Segalini (2020)','Location','southwest')

sgtitle('Rotor averaged wind speed evolution of the first turbine');




