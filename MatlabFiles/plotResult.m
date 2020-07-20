clear
clc
close all

%% Parameters
m = 2000;
g = 9.81;
crr = 0.006;
CdA = 0.6;
ds = 10;
eta_trans = 0.98;
eta_dc = 0.99;
alpha0 = 785.0;
alpha1 = -43e-4;
alpha2 = 5.6e-7;
beta0 = 41e-7;
penalty = 11e4;

Vmax = 200/3.6;
Vmin = 0.0;
Fmax = 3e3;
Fmin = -3e3;
PAmax = 6e4;
PDmax = -6e4;

angle = 0;

%% Scenario 1
% Problem sizes 
Nx = 200;
Nu = 500;
Nhrz = 200;

% The result of using static grid without calibration
Vo_nocalibrate = importdata('./SpeedSolverData/Scenario_1/Vo_nocalibrate_1.mat');
Fo_nocalibrate = importdata('./SpeedSolverData/Scenario_1/Fo_nocalibrate_1.mat');
MinCost_nocalibrate = importdata('./SpeedSolverData/Scenario_1/MinCost_nocalibrate_1.mat');

% The result of using static grid with calibration
Vo_static = importdata('./SpeedSolverData/Scenario_1/Vo_static_1.mat');
Fo_static = importdata('./SpeedSolverData/Scenario_1/Fo_static_1.mat');
MinCost_static = importdata('./SpeedSolverData/Scenario_1/MinCost_static_1.mat');

% The result of using dynamic grid
Vo_dynamic = importdata('./SpeedSolverData/Scenario_1/Vo_dynamic_1.mat');
Fo_dynamic = importdata('./SpeedSolverData/Scenario_1/Fo_dynamic_1.mat');
MinCost_dynamic = importdata('./SpeedSolverData/Scenario_1/MinCost_dynamic_1.mat');

% Boundary Info
upperBound = importdata('./SpeedSolverData/Scenario_1/upperBound_1.mat');
lowerBound = importdata('./SpeedSolverData/Scenario_1/lowerBound_1.mat');
upperActual = importdata('./SpeedSolverData/Scenario_1/upperActual_1.mat');
lowerActual = importdata('./SpeedSolverData/Scenario_1/lowerActual_1.mat');

% Environmental Info
envFactor = struct;
envFactor.Vmax_env = zeros((Nhrz+1), 1);
envFactor.Vmin_env = zeros((Nhrz+1), 1);
envFactor.Angle_env = zeros((Nhrz+1), 1);

Vmax_GPS = [130/3.6 80/3.6 50/3.6 80/3.6];
Vmin_GPS = [60/3.6 30/3.6 10/3.6 30/3.6];
Angle_GPS = [0.0 0.0 0.0 0.0];
changePoint = [30 100 150];

envFactor.endBlock = [changePoint Nhrz];

for i = 1:(Nhrz+1)
    if i<=changePoint(1)
        envFactor.Vmax_env(i) = Vmax_GPS(1);
        envFactor.Vmin_env(i) = Vmin_GPS(1);
        envFactor.Angle_env(i) = Angle_GPS(1);
    elseif i<=changePoint(2)
        envFactor.Vmax_env(i) = Vmax_GPS(2);
        envFactor.Vmin_env(i) = Vmin_GPS(2);
        envFactor.Angle_env(i) = Angle_GPS(2);
    elseif i<=changePoint(3)
        envFactor.Vmax_env(i) = Vmax_GPS(3);
        envFactor.Vmin_env(i) = Vmin_GPS(3);
        envFactor.Angle_env(i) = Angle_GPS(3);
    else
        envFactor.Vmax_env(i) = Vmax_GPS(4);
        envFactor.Vmin_env(i) = Vmin_GPS(4);
        envFactor.Angle_env(i) = Angle_GPS(4);
    end 
end

%% Scenario 2
% Problem sizes 
Nx = 200;
Nu = 500;
Nhrz = 100;

% The result of using static grid without calibration
Vo_nocalibrate = importdata('./SpeedSolverData/Scenario_2/Vo_nocalibrate_2.mat');
Fo_nocalibrate = importdata('./SpeedSolverData/Scenario_2/Fo_nocalibrate_2.mat');
MinCost_nocalibrate = importdata('./SpeedSolverData/Scenario_2/MinCost_nocalibrate_2.mat');

% The result of using static grid with calibration
Vo_static = importdata('./SpeedSolverData/Scenario_2/Vo_static_2.mat');
Fo_static = importdata('./SpeedSolverData/Scenario_2/Fo_static_2.mat');
MinCost_static = importdata('./SpeedSolverData/Scenario_2/MinCost_static_2.mat');

% The result of using dynamic grid
Vo_dynamic = importdata('./SpeedSolverData/Scenario_2/Vo_dynamic_2.mat');
Fo_dynamic = importdata('./SpeedSolverData/Scenario_2/Fo_dynamic_2.mat');
MinCost_dynamic = importdata('./SpeedSolverData/Scenario_2/MinCost_dynamic_2.mat');

% Boundary Info
upperBound = importdata('./SpeedSolverData/Scenario_2/upperBound_2.mat');
lowerBound = importdata('./SpeedSolverData/Scenario_2/lowerBound_2.mat');
upperActual = importdata('./SpeedSolverData/Scenario_2/upperActual_2.mat');
lowerActual = importdata('./SpeedSolverData/Scenario_2/lowerActual_2.mat');

% Environmental Info
envFactor = struct;
envFactor.Vmax_env = zeros((Nhrz+1), 1);
envFactor.Vmin_env = zeros((Nhrz+1), 1);
envFactor.Angle_env = zeros((Nhrz+1), 1);

Vmax_GPS = [80/3.6];
Vmin_GPS = [0/3.6];
Angle_GPS = [0.0];
changePoint = [];

envFactor.endBlock = [changePoint Nhrz];

for i = 1:(Nhrz+1)
        envFactor.Vmax_env(i) = Vmax_GPS(1);
        envFactor.Vmin_env(i) = Vmin_GPS(1);
        envFactor.Angle_env(i) = Angle_GPS(1);
end

%% Set up Upper and Lower bounds of speed
[a, N] = size(envFactor.endBlock);
[b, M] = size(Vmax_GPS);

endBlock = envFactor.endBlock;
endBlock(end) = endBlock(end) + 1;
turningPoint = [1 endBlock];

%% Verify the speed trajectory (based on the given control policy and inital speed)
% Parameter Settings (local copy)
% ActualSpeedTrajectory = zeros(1, Nhrz+1);
% ActualSpeedTrajectory(1) = Vo_nocalibrate(1);
% 
% for i = 1:Nhrz
%     Xk_plus_1 = sqrt((2*ds/m)*Fo_nocalibrate(i) + (1 - 2*ds*CdA/m)*(ActualSpeedTrajectory(i))^2 - 2*ds*g*(sin(angle)+crr*cos(angle)));
%     ActualSpeedTrajectory(i+1) = Xk_plus_1;
% end

%% Plot
% Generate State Vector
delta = (Vmax - Vmin) / (Nx - 1);
StateVector = zeros(1, Nx);

for i = 2:Nx
    StateVector(i) = StateVector(i - 1) + delta;
end

hold on
% 
for i = 0:10
    plot(i*ds, (StateVector')*3.6, "k.",'MarkerSize',8);
end

% Plot speed trajectory, boundary, etc.
grid on;
%plot((0:Nhrz)*ds, ActualSpeedTrajectory(1:NHrz+1)*3.6,'y','LineWidth',1.2);
line(1) = plot((0:Nhrz)*ds, Vo_nocalibrate(1:Nhrz+1)*3.6,'-','LineWidth',1.2, 'Color', [0.9290, 0.6940, 0.1250]);
line(2) = plot((0:Nhrz)*ds, Vo_static(1:Nhrz+1)*3.6,'-','LineWidth',1.2, 'Color', [0.8500, 0.3250, 0.0980]);
line(3) = plot((0:Nhrz)*ds, Vo_dynamic(1:Nhrz+1)*3.6,'-','LineWidth',1.2, 'Color', [0, 0.4470, 0.7410]);

for i = 1:N
    yUpper = [Vmax_GPS(i), Vmax_GPS(i)];
    yLower = [Vmin_GPS(i), Vmin_GPS(i)];
    xIndex = [turningPoint(i), (turningPoint(i+1) - 1)];
    line(4) = plot(xIndex*ds, yUpper*3.6, 'LineWidth',2, 'Color', [0.25, 0.25, 0.25]);
    line(5) = plot(xIndex*ds, yLower*3.6, 'LineWidth',2, 'Color', [0.25, 0.25, 0.25]);
end
line(6) = plot((0:Nhrz)*ds, upperBound(1:Nhrz+1)*3.6,'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
line(7) = plot((0:Nhrz)*ds, lowerBound(1:Nhrz+1)*3.6,'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
line(8) = plot((0:Nhrz)*ds, upperActual(1:Nhrz+1)*3.6,'r.','MarkerSize',5);
line(9) = plot((0:Nhrz)*ds, lowerActual(1:Nhrz+1)*3.6,'r.','MarkerSize',5);

xlabel('Distance (m)');
ylabel('Speed (km/h)');
legend(line([1 2 3 4 6 8]),  'Speed trajectpry without boundary calibration', 'Speed trajectpry with boundary calibration', 'Speed trajectpry with using nonlinear state grid', 'Legal speed limits', 'Boundary line', 'Calibration')
hold off;

disp('The Minimum cost of static without calibration:')
disp(MinCost_nocalibrate)
disp('The Minimum cost of static with calibration:')
disp(MinCost_static)
disp('The Minimum cost of dynamic:')
disp(MinCost_dynamic)