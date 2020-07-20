%% Mode settings
boundaryFlag = 0;
gridFlag = 0;
checkTrajectoryFlag = 0;

%% Extract the results from Simulink
Fo_copy = Fo(1,:);
Vo_copy = Vo(1,:);
Vo_copy = [V0 Vo_copy];
Qo_copy = Qo(1,:);
To_copy = To(1,:);
To_copy = [T0 To_copy];
MinCost = MinCost(1);

%% Boundary
if (boundaryFlag)
    upperSpeedBound_copy = upperSpeedBound(1,:);
    lowerSpeedBound_copy = lowerSpeedBound(1,:);
    upperSpeedActual_copy = [X0 upperSpeedActual(1,:)];
    lowerSpeedActual_copy = [X0 lowerSpeedActual(1,:)];
    upperTempBound_copy = upperTempBound(1,:);
    lowerTempBound_copy = lowerTempBound(1,:);
end


%% Local Copy of grid sizes from 'createSfunc.m'
NHrz = solverinput.GridSize.Nhrz;
Nhrz = solverinput.GridSize.Nhrz;
Nv = solverinput.GridSize.Nv;
ResThermal = solverinput.GridSize.ResThermal;
Nt = solverinput.GridSize.Nt;

%% Set up Upper and Lower bounds of speed
[a, N] = size(envFactor.endBlock);
[b, M] = size(Vmax_GPS);

endBlock = envFactor.endBlock;
endBlock(end) = endBlock(end) + 1;
turningPoint = [1 endBlock];

%% Verify the speed trajectory (based on the given control policy and inital speed)
% Parameter Settings (local copy)
m = modelPara.m;
g = modelPara.g;
CdA = modelPara.CdA;
crr = modelPara.crr;
ds = modelPara.ds;
eta_trans = modelPara.eta_trans;
eta_dc = modelPara.eta_dc;
alpha0 = modelPara.alpha0;
alpha1 = modelPara.alpha1;
alpha2 = modelPara.alpha2;
beta0 = modelPara.beta0;
penalty = modelPara.speedPenalty;

angle = 0.0;

ActualSpeedTrajectory = zeros(1, NHrz+1);
ActualSpeedTrajectory(1) = V0;

for i = 1:NHrz
    Xk_plus_1 = sqrt((2*ds/m)*Fo_copy(i) + (1 - 2*ds*CdA/m)*(ActualSpeedTrajectory(i))^2 - 2*ds*g*(sin(angle)+crr*cos(angle)));
    ActualSpeedTrajectory(i+1) = Xk_plus_1;
end

%% Verify the Total cost (based on the given speed trajectory and control policy)
totalCost = 0;

for i = 1:NHrz
    
    P_wh = Vo_copy(i)*Fo_copy(i);
    
    if Fo_copy(i)>0
        % Acceleration
        P_m = P_wh/eta_trans;
        P_inv = ((1-alpha1)-sqrt((alpha1-1)^2 - 4*alpha2*(alpha0+P_m)))/(2*alpha2);
        P_dc = P_inv/eta_dc;
        P_batt = (1 - sqrt(1-4*beta0*P_dc))/(2*beta0);
        disp(i)
        disp(P_batt)
    else
        % Deceleration
        P_m = P_wh*eta_trans;
        P_inv = ((1-alpha1)-sqrt((alpha1-1)^2 - 4*alpha2*(alpha0+P_m)))/(2*alpha2);
        P_dc = P_inv*eta_dc;
        P_batt = (1 - sqrt(1-4*beta0*P_dc))/(2*beta0);
    end
    
    dt = 2*ds/(Vo_copy(i+1) + Vo_copy(i));
    
    totalCost = totalCost + (P_batt + penalty)*dt;
end

disp('Actual Cost based on the given trajectories:')
disp(totalCost)
disp('The Minimum cost is:')
disp(MinCost)

%% Plot - Speed
% Generate State Vector - Speed
delta = (solverinput.Constraint.Vmax - solverinput.Constraint.Vmin) / (Nv - 1);
SpeedVector = zeros(1, Nv);

for i = 2:Nv
    SpeedVector(i) = SpeedVector(i - 1) + delta;
end

figure(1)
hold on

if (gridFlag)
    for i = 0:10
        plot(i*ds, (SpeedVector')*3.6, "k.",'MarkerSize',8);
    end
end

% Plot speed trajectory, boundary, etc.
grid on;
line(1) = plot((0:Nhrz)*ds, Vo_copy(1:Nhrz+1)*3.6,'-','LineWidth',1.2, 'Color', [0.8500, 0.3250, 0.0980]);


for i = 1:N
    yUpper = [Vmax_GPS(i), Vmax_GPS(i)];
    yLower = [Vmin_GPS(i), Vmin_GPS(i)];
    xIndex = [turningPoint(i), (turningPoint(i+1) - 1)];
    line(2) = plot(xIndex*ds, yUpper*3.6, 'LineWidth',2, 'Color', [0.25, 0.25, 0.25]);
    line(3) = plot(xIndex*ds, yLower*3.6, 'LineWidth',2, 'Color', [0.25, 0.25, 0.25]);
end


if (checkTrajectoryFlag)
    line(4) = plot((0:Nhrz)*ds, ActualSpeedTrajectory(1:NHrz+1)*3.6,'LineWidth',1.2, 'Color', [0.9290, 0.6940, 0.1250]);
end

if(boundaryFlag)
    line(5) = plot((0:Nhrz)*ds, upperSpeedBound_copy(1:Nhrz+1)*3.6,'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
    line(6) = plot((0:Nhrz)*ds, lowerSpeedBound_copy(1:Nhrz+1)*3.6,'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
    line(7) = plot((0:Nhrz)*ds, upperSpeedActual_copy(1:Nhrz+1)*3.6,'r.','MarkerSize',5);
    line(8) = plot((0:Nhrz)*ds, lowerSpeedActual_copy(1:Nhrz+1)*3.6,'r.','MarkerSize',5);
end

xlabel('Distance (m)');
ylabel('Speed (km/h)');

%legend(line([1 2 4 5 7]),  'Calculated speed trajectpry', 'Legal speed limits', 'Actual speed trajectory','Calculated Boundary line', 'Calibration')
legend(line([1 2]),  'Calculated speed trajectpry', 'Legal speed limits')

hold off;

%% Plot - Thermal
% Generate State Vector - Thermal
delta = (solverinput.Constraint.Tmax - solverinput.Constraint.Tmin) / (Nt - 1);
ThermalVector = zeros(1, Nt);
ThermalVector(1) = solverinput.Constraint.Tmin;

for i = 2:Nt
    ThermalVector(i) = ThermalVector(i - 1) + delta;
end

figure(2)
hold on

if (gridFlag)
    for i = 0:10
        plot(i*ds, (ThermalVector')*3.6, "r.",'MarkerSize',8);
    end
end


% Plot thermal trajectory, boundary, etc.
grid on;
line(1) = plot((0:ResThermal)*10*ds, To_copy(1:ResThermal+1),'-','LineWidth',1.2, 'Color', [0.8500, 0.3250, 0.0980]);

for i = 1:N
    yRequired = [T_required(i), T_required(i)];
    xIndex = [turningPoint(i), (turningPoint(i+1) - 1)];
    line(2) = plot(xIndex*ds, yRequired, 'LineWidth',2, 'Color', [0.25, 0.25, 0.25]);
end

if(boundaryFlag)
    line(5) = plot((0:ResThermal)*10*ds, upperTempBound_copy(1:ResThermal+1),'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
    line(6) = plot((0:ResThermal)*10*ds, lowerTempBound_copy(1:ResThermal+1),'LineWidth',1.2, 'Color', [0.4660, 0.6740, 0.1880]);
    line(7) = plot((0:Nhrz)*ds, upperTempActual_copy(1:Nhrz+1)*3.6,'r.','MarkerSize',5);
    line(8) = plot((0:Nhrz)*ds, lowerTempActual_copy(1:Nhrz+1)*3.6,'r.','MarkerSize',5);
end

xlabel('Distance (m)');
ylabel('Cabin Temperature (celsius)');

legend(line([1 2]),  'Calculated thermal trajectpry', 'Required temperature')

hold off;

