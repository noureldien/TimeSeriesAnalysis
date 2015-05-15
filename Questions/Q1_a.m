%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1-a: Auto Regression
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
arOrder = 3;
[yEstm, arVariance, arResidual, arParams] = autoRegression(yTrain, arOrder);

% % plot the AR estimates vs the training
% figure(1); clf;
% subplot(2,1,1);
% hold on;
% grid on;
% box on;
% axis([0 120 -inf inf]);
% plot(yTrain, 'b', 'LineWidth', 1);
% plot(yEstm, 'r', 'LineWidth', 1);
% xlabel('Time (month)', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Index Prediction using AR(3) Model', 'FontSize', 16);
% plot_legend = legend('Actual', 'Predicted', 'Location', 'SE');
% set(plot_legend, 'FontSize', 10);
% subplot(2,1,2);
% hold on;
% grid on;
% box on;
% axis([0 120 -inf inf]);
% plot(abs(arResidual), 'b', 'LineWidth', 1);
% xlabel('Time (month)', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Prediction Error (Absolute)', 'FontSize', 16);








