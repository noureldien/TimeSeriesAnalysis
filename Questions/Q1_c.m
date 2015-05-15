%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1-a: Lasso Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

% load data
load('data/dataFeatures');
load('data/yTrain');
load('data/yEstmKf');

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

% lassoTarget is the residue of the Kalman filter
% make sure to normalize it
lassoTarget = yTrain - yEstmKf;
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

% plot the errors of lasso
figure(1); clf;
hold on;
grid on;
box on;
plot(lassoLambda, lassoErrors, 'LineWidth', 2);
xlabel('Lasso Regulariser', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
title('Error (Absolute) of Lasso Regression', 'FontSize', 16);

% use the weights of the smallest error
% it is chosen imperically
lassoResult = normFeatures * lassoWeights(:, 81);

% plot the weights
figure(2); clf;
colorMap = lines(10);
subplot(2,1,1);
hold on;
grid on;
box on;
plot(lassoLambda, lassoWeights', 'LineWidth', 2);
ylabel('Weight Value', 'FontSize', 16);
title('Decay of Feature Weights in Lasso Regression', 'FontSize', 16);
plot_legend = legend('BER','OIL','PMI','INCOME','PROFIT','POP','UNEMP', 'Location', 'SE');
set(plot_legend, 'FontSize', 10);
subplot(2,1,2);
colorMap = lines(10);
hold on;
grid on;
box on;
plot(lassoLambda, lassoWeights', 'LineWidth', 2);
xlabel('Lasso Regulariser', 'FontSize', 16);
ylabel('Weight Value', 'FontSize', 16);
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







