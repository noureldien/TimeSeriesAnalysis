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

% plot the AR estimates vs the training
coloMap = lines(30);
colorGreen = [0 0.7 0.2];
figure(1); clf;
subplot(2,1,1);
hold on;
grid on;
box on;
plot1 = plot(yEstmAR, 'LineWidth', 1, 'Color', coloMap(1,:));
plot2 = plot(yEstmKf, 'LineWidth', 1, 'Color', coloMap(2,:));
plot3 = plot(yTrain, 'LineWidth', 1, 'Color', colorGreen);
xlabel('Time (month)', 'FontSize', 16);
ylabel('Stock Value', 'FontSize', 16);
title('Index Prediction', 'FontSize', 16);
plot_legend = legend([plot3, plot1, plot2], {'Actual', 'AR Prediction', 'KF Prediction'}, 'Location', 'SE');
set(plot_legend, 'FontSize', 12);
subplot(2,1,2);
hold on;
grid on;
box on;
plot(errorAR, 'LineWidth', 1, 'Color', coloMap(1,:));
plot(errorKm, 'LineWidth', 1, 'Color', coloMap(2,:));
xlabel('Time (month)', 'FontSize', 16);
ylabel('Stock Value', 'FontSize', 16);
title('Prediction Error', 'FontSize', 16);
plot_legend = legend('AR Prediction', 'KF Prediction', 'Location', 'SE');
set(plot_legend, 'FontSize', 12);









