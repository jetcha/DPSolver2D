%% Results
computeTime1 = [];
computeTime2 = [];
computeTime3 = [];
computeTime4 = [];
computeTime5 = [];

totalCost1 = [4.950e+07 4.913e+07 4.908e+07 4.906e+07 4.906e+07 4.906e+07 4.906e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07 4.905e+07];
          
totalCost2 = [4.950e+07 4.824e+07 4.743e+07 4.718e+07 4.709e+07 4.706e+07 4.706e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07 4.705e+07];
          
totalCost3 = [4.950e+07 4.784e+07 4.690e+07 4.660e+07 4.649e+07 4.641e+07 4.639e+07 4.637e+07 4.637e+07 4.637e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07 4.636e+07];
          
totalCost4 = [4.950e+07 4.813e+07 4.700e+07 4.676e+07 4.653e+07 4.644e+07 4.635e+07 4.631e+07 4.628e+07 4.627e+07 4.626e+07 4.625e+07 4.624e+07 4.624e+07 4.624e+07 4.624e+07 4.624e+07 4.624e+07 4.624e+07 4.624e+07];
          
totalCost5 = [4.950e+07 4.815e+07 4.809e+07 4.757e+07 4.731e+07 4.708e+07 4.697e+07 4.660e+07 4.647e+07 4.643e+07 4.638e+07 4.637e+07 4.636e+07 4.633e+07 4.630e+07 4.630e+07 4.629e+07 4.629e+07 4.628e+07 4.628e+07];

iteration = 1:20;

%% Results 2
computeTime1 = [];
computeTime2 = [];
computeTime3 = [];
computeTime4 = [];
computeTime5 = [];

totalCost1 = [4.458e+07 4.446e+07 4.441e+07 4.434e+07 4.427e+07 4.427e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07 4.426e+07];  
          
totalCost2 = [4.458e+07 4.228e+07 4.208e+07 4.175e+07 4.167e+07 4.164e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07 4.163e+07];  
          
totalCost3 = [4.458e+07 4.240e+07 4.138e+07 4.104e+07 4.087e+07 4.080e+07 4.078e+07 4.076e+07 4.076e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07 4.075e+07];      
          
totalCost4 = [4.458e+07 4.205e+07 4.124e+07 4.070e+07 4.047e+07 4.036e+07 4.024e+07 4.020e+07 4.017e+07 4.014e+07 4.011e+07 4.010e+07 4.009e+07 4.008e+07 4.008e+07 4.007e+07 4.007e+07 4.007e+07 4.007e+07 4.007e+07];
          
totalCost5 = [4.458e+07 4.239e+07 4.184e+07 4.139e+07 4.115e+07 4.083e+07 4.072e+07 4.057e+07 4.050e+07 4.049e+07 4.044e+07 4.038e+07 4.034e+07 4.032e+07 4.030e+07 4.029e+07 4.025e+07 4.024e+07 4.023e+07 4.021e+07];

iteration = 1:20;

totalCost1 = [4.920e+07 4.758e+07 4.681e+07 4.651e+07 4.631e+07 4.621e+07 4.618e+07 4.616e+07 4.615e+07 4.615e+07];

totalCost1 = [4.279e+07 4.1571e+07 4.089e+07 4.059e+07 4.044e+07 4.036e+07 4.031e+07 4.027e+07 4.025e+07 4.024e+07];  
%% Color
colorArray1 = [0, 0.4470, 0.7410];
colorArray2 = [0.8500, 0.3250, 0.0980];
colorArray3 = [0.4660, 0.6740, 0.1880];
colorArray4 = [0.9290, 0.6940, 0.1250];
colorArray5 = [0.4940, 0.1840, 0.5560];

%% Plot

figure(3)
hold on 
grid on;
line(1) = plot(iteration, totalCost1,'-o','LineWidth',1.2, 'Color', colorArray1);
% line(2) = plot(iteration, computeTime2,'-o','LineWidth',1.2, 'Color', colorArray2);
% line(3) = plot(iteration, computeTime3,'-o','LineWidth',1.2, 'Color', colorArray3);
% line(4) = plot(iteration, computeTime4,'-o','LineWidth',1.2, 'Color', colorArray4);
% line(5) = plot(iteration, computeTime5,'-o','LineWidth',1.2, 'Color', colorArray5);
title('Computational time with different N_v and N_f')
xlabel('Iteration');
ylabel('Computational time (s)');
%ylim([24 27])
%xlim([1 20])
%legend(line([1 2 3 4 5]),  {'\Delta s = 10', '\Delta s = 20', '\Delta s = 30', '\Delta s = 40', '\Delta s = 50'},'Location','northwest')
%legend(line([1 2 3 4 5]),  {'N_v = N_f = 11', 'N_v = N_f = 21', 'N_v = N_f = 31', 'N_v = N_f = 41', 'N_v = N_f = 51'},'Location','northwest')
hold off;

figure(4)
hold on 
grid on;
line(1) = plot(iteration, totalCost1,'-o','LineWidth',1.2, 'Color', colorArray1);
%line(2) = plot(iteration, totalCost1,'-o','LineWidth',1.2, 'Color', colorArray2);
% line(3) = plot(iteration, totalCost3,'-o','LineWidth',1.2, 'Color', colorArray3);
% line(4) = plot(iteration, totalCost4,'-o','LineWidth',1.2, 'Color', colorArray4);
% line(5) = plot(iteration, totalCost5,'-o','LineWidth',1.2, 'Color', colorArray5);
title('Total cost over iterations')
xlabel('Iteration');
ylabel('Total cost');
%ylim([24 27])
%xlim([1 20])
%legend(line([1 2 3 4 5]),  {'\Delta s = 10', '\Delta s = 20', '\Delta s = 30', '\Delta s = 40', '\Delta s = 50'},'Location','northeast')
%legend(line([1 2 3 4 5]),  {'\lambda = 0.25', '\lambda = 0.33', '\lambda = 0.5', '\lambda = 0.67', '\lambda = 0.83'},'Location','northeast')
hold off;

% 
% figure(5)
% hold on 
% grid on;
% line(1) = plot(iteration, speedIDAE,'-o','LineWidth',1.2, 'Color', colorArray);
% title('IDAE(speed) over iterations')
% xlabel('Iteration');
% ylabel('IDAE - speed');
% %ylim([24 27])
% %xlim([0 Nhrz*ds])
% %legend(line([1 2 3]),  {'Solver Result', 'Model based result', 'Required temperature'},'Location','northwest')
% hold off;
% 
% figure(6)
% hold on 
% grid on;
% line(1) = plot(iteration, tempIDAE,'-o','LineWidth',1.2, 'Color', colorArray);
% title('IDAE(temperature) over iterations')
% xlabel('Iteration');
% ylabel('IDAE - temperature');
% %ylim([24 27])
% %xlim([0 Nhrz*ds])
% %legend(line([1 2 3]),  {'Solver Result', 'Model based result', 'Required temperature'},'Location','northwest')
% hold off;