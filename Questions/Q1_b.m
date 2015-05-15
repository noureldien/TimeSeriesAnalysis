%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1-a: Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

% load data
load('data/stockIndex.mat');

% we depend on the 4th column, which is the close price
y = stockIndex(:,4);

% get stock index (monthly) from 1/May/1997 to 1/May/2007
% and this will be the training

% for testing, we take from 2/May/2007 to 1/May/2015
ySplitIndex = 121;
yTrain = y(1:ySplitIndex);
yTest = y(ySplitIndex+1:end);

% get the V, which is the observation noise
% using the Auto-Regression model
kfOrder = 3;
kfWindow = 3;
[yEstmAR, kfS] = autoRegression(yTrain, kfOrder);

% get estimate using Kalman filter
yEstmKf = kalmanFilter(yTrain, kfOrder, kfS, kfWindow);

errorAR = abs(yTrain - yEstmAR);
errorKm = abs(yTrain - yEstmKf);

% % plot the AR/Kalman estimates vs the training
% coloMap = lines(30);
% colorGreen = [0 0.7 0.2];
% figure(2); clf;
% subplot(2,1,1);
% hold on;
% grid on;
% box on;
% axis([0 120 700 1600]);
% plot1 = plot(yEstmAR, 'LineWidth', 1, 'Color', 'r');
% plot2 = plot(yEstmKf, 'LineWidth', 1, 'Color', 'b');
% plot3 = plot(yTrain, 'LineWidth', 1, 'Color', colorGreen);
% xlabel('Time (month)', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Index Prediction using Kalman and AR', 'FontSize', 16);
% plot_legend = legend([plot3, plot1, plot2], {'Actual', 'AR Prediction', 'KF Prediction'}, 'Location', 'SE');
% set(plot_legend, 'FontSize', 10);
% subplot(2,1,2);
% hold on;
% grid on;
% box on;
% axis([0 120 -inf 250]);
% plot(abs(errorAR), 'LineWidth', 1, 'Color', 'r');
% plot(abs(errorKm), 'LineWidth', 1, 'Color', 'b');
% xlabel('Time (month)', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Prediction Error (Absolute)', 'FontSize', 16);
% plot_legend = legend('AR Prediction', 'KF Prediction', 'Location', 'NE');
% set(plot_legend, 'FontSize', 10);

% figure(3); clf;
% hold on;
% grid on;
% box on;
% boxplot([abs(errorAR), abs(errorKm)]);
% ylim([-20 200]);
% xlabel('AR             KF', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Error (Absolute)', 'FontSize', 16);

% % plot effect of small and high QR ratio
% colorGreen = [0 0.7 0.2];
% figure(2); clf;
% hold on;
% grid on;
% box on;
% axis([0 120 700 1600]);
% plot1 = plot(yEstmKf1, 'LineWidth', 1, 'Color', 'r');
% plot2 = plot(yEstmKf2, 'LineWidth', 1, 'Color', 'b');
% plot3 = plot(yTrain, 'LineWidth', 1, 'Color', colorGreen);
% xlabel('Time (month)', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Index Prediction using Kalman and AR', 'FontSize', 16);
% plot_legend = legend([plot3, plot1, plot2], {'Actual', 'Avg. RQ Ratio', 'Big RQ Ratio'}, 'Location', 'SW');
% set(plot_legend, 'FontSize', 10);











