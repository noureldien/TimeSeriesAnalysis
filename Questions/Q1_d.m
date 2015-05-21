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

% normalize the regression target
yTrain = yTrain - mean(yTrain);
yTrain = yTrain / std(yTrain);

yEstmAR = yEstmAR - mean(yEstmAR);
yEstmAR = yEstmAR / std(yEstmAR);

yEstmKf = yEstmKf - mean(yEstmKf);
yEstmKf = yEstmKf / std(yEstmKf);

% estimate of the index using LagLasso
yEstmLL = yTrain;
lagWindow = 10;

lambdas = [];

%lassoTarget is the sliding window of stock index
for i=lagWindow+1:size(yTrain, 1)
    i
    idxS = i-lagWindow;
    idxE = i-1;
    lassoTarget = yTrain(idxS:idxE);
    lassoFeatures = normFeatures(idxS:idxE,:);
       
    % lasso regression
    [lassoWeights, lassoInfo] = lasso(lassoFeatures, lassoTarget);
    lassoLambda = lassoInfo.Lambda;
    
    % find the best set of weights by calculating the error
    lassoErrors = zeros(length(lassoInfo.Lambda),1);
    for j=1:length(lassoErrors)
        lassoResult = lassoFeatures * lassoWeights(:, j);
        lassoErrors(j) = mean(abs(lassoResult-lassoTarget));
    end
    
    [~,idx] = sort(lassoErrors);
    lambdas = [lambdas, lassoLambda(idx(1))];
    % estimate current stock index
    yEstmLL(i) = normFeatures(i,:) * lassoWeights(:, idx(1));
end

% plot regression error
figure(1); clf;
hold on;
grid on;
box on;
plot(abs(yEstmLL-yTrain), 'LineWidth', 2);
xlabel('Lasso Regulariser', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
title('Error (Absolute) of Lasso Regression', 'FontSize', 16);

% plot the AR/Kalman/LagLasso estimates vs the training
coloMap = lines(30);
colorGreen = [0 0.7 0.2];
figure(3); clf;
hold on;
grid on;
box on;
plot(yTrain, 'LineWidth', 1, 'Color', 'k');
plot(yEstmAR, 'LineWidth', 1, 'Color', 'r');
plot(yEstmKf, 'LineWidth', 1, 'Color', 'b');
plot(yEstmLL, 'LineWidth', 1, 'Color', colorGreen);
xlabel('Time (month)', 'FontSize', 16);
ylabel('Value', 'FontSize', 16);
title('Index Prediction using Kalman and AR', 'FontSize', 16);