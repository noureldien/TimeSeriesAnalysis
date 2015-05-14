function [yEstm, arVariance, arResidual, arParams] = autoRegression(yTrain, arOrder)

% calcuate the V, which is the observation noise, as the following:
% build the auto-regression model (2nd/3rd) and calcuate the residual
% the variance of the residual can be used as the observation noise

arOpt = arOptions('Approach', 'ls', 'EstCovar', true, 'Window', 'pow');
[arModel, arParams] = ar(yTrain,arOrder,arOpt);
yEstm = compare(yTrain,arModel,1);
arResidual = yTrain - yEstm;

% variance of resudual is the obersvation noise of kalman
% can be optained by other method:
% arVariance = arModel.NoiseVariance;
% arVariance = var(arResidual);
arVariance = std(arResidual);

end