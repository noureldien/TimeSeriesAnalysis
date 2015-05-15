%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1-a: LagLasso Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

% load data
load('data/dataFeatures');
load('data/yTrain');
load('data/yEstmKf');
load('data/yEstmAR');

% normalize the predictors, i.e the features
normFeatures = zeros(size(dataFeatures));
for i=1:size(normFeatures,2)
    normFeatures(:,i) = dataFeatures(:,i) - mean(dataFeatures(:,i));
    normFeatures(:,i) = normFeatures(:,i) / std(normFeatures(:,i));
end

% build a lasso regression model to see which of these
% features is of a great important, i.e. which of these
% micro-ecoonomic predictor has a greate influence
% in predicting the stock index

% lassoTarget is the sliding window of stock index
y
% for i=window+1:size(yTrain, 1)
%     %lassoWeights = lassoRegression(dataFeatures, yTrain, lassoTaw);
%     idxS = i-30;
%     idxE = i-1;
%     target = yTrain(idxS:idxE);
%     features = normFeatures(idxS:idxE,:);
%     %[lassoWeights, lassoInfo] = lasso(features, target, ...
%     %                               'Alpha', lassoTaw, 'Lambda', ones(window, 1));
%     lassoWeights = lassoRegression(features, target, lassoTaw);
%     lassoWeights = lassoRegression(features.*lassoWeights', target, lassoTaw);
%     result = features*lassoWeights;
% end


lassoTarget = yTrain;

lassoTarget = lassoTarget - mean(lassoTarget);
lassoTarget = lassoTarget / std(lassoTarget);

% lasso regression
[lassoWeights,lassoInfo] = lasso(normFeatures,lassoTarget);
lassoLambda = lassoInfo.Lambda;

% find the best set of weights by calculating the error
lassoErrors = zeros(length(lassoInfo.Lambda),1);
for i=1:length(lassoErrors)
lassoResult = normFeatures * lassoWeights(:, i);
lassoErrors(i) = mean(abs(lassoResult-lassoTarget));
end

[~,idx] = sort(lassoErrors);
% use the weights of the smallest error
% it is chosen imperically
lassoResult = normFeatures * lassoWeights(:, idx(1));

% plot the errors of lasso
figure(1); clf;
hold on;
grid on;
box on;
plot(lassoLambda, lassoErrors, 'LineWidth', 2);
xlabel('Lasso Regulariser', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
title('Error (Absolute) of Lasso Regression', 'FontSize', 16);

% plot the weights
figure(2); clf;
colorMap = lines(10);
hold on;
grid on;
box on;
plot(lassoLambda, lassoWeights', 'LineWidth', 2);
xlabel('Lasso Regulariser', 'FontSize', 16);
ylabel('Weight Value', 'FontSize', 16);
title('Decay of Feature Weights in Lasso Regression', 'FontSize', 16);
plot_legend = legend('BER','OIL','PMI','INCOME','PROFIT','POP','UNEMP', 'Location', 'SE');
set(plot_legend, 'FontSize', 10);

return;

% plot the AR/Kalman estimates vs the training
coloMap = lines(30);
colorGreen = [0 0.7 0.2];
figure(3); clf;
hold on;
grid on;
box on;
plot1 = plot(lassoTarget, 'LineWidth', 1, 'Color', 'b');
plot2 = plot(lassoResult, 'LineWidth', 1, 'Color', 'r');
xlabel('Time (month)', 'FontSize', 16);
ylabel('Value', 'FontSize', 16);
title('Index Prediction using Kalman and AR', 'FontSize', 16);

return;

% estimates of the stock index after sparse regression
yEstmLassoTrain = normFeatures*lassoWeights;
%yEstmLassoTrain = yEstmLassoTrain+1200;

% plot the AR/Kalman estimates vs the training
coloMap = lines(30);
colorGreen = [0 0.7 0.2];
figure(2); clf;
%subplot(2,1,1);
hold on;
grid on;
box on;
axis([0 120 700 1600]);
plot1 = plot(yEstmAR, 'LineWidth', 1, 'Color', 'r');
plot2 = plot(yEstmKf, 'LineWidth', 1, 'Color', 'b');
plot3 = plot(yTrain, 'LineWidth', 1, 'Color', colorGreen);
plot4 = plot(yEstmLassoTrain, 'LineWidth', 1, 'Color', 'k');
xlabel('Time (month)', 'FontSize', 16);
ylabel('Value', 'FontSize', 16);
title('Index Prediction using Kalman and AR', 'FontSize', 16);
%plot_legend = legend([plot3, plot1, plot2], {'Actual', 'AR Prediction', 'KF Prediction'}, 'Location', 'SE');
%set(plot_legend, 'FontSize', 10);
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











