function [yEstm] = kalmanFilter(yTrain, order, s, window)

% Generate synthetic and ideal data and fit a DLM smoother.
% The command |dlmgensys| generates system evolution matrix |G|
% and obs operator |F|. We use local level and trend model with 12 seasons.
% ns = number of seasons
sysOptions = struct('order',1,'fullseas',1,'ns',window);
[G,F] = dlmgensys(sysOptions);
% p = number of data sets = 1,
% m = number of internal states = 13
[p,m] = size(F);

% Generate data.
n = size(yTrain,1);

V = ones(n,p)*s;

% initial state must be drawn from ~N(0, Vo)
x0 = zeros(window+1, 1);
for i=window+1
    x0(i) = normrnd(0,s);
end

% initial state uncertainty
C0 = eye(m)*s^2;

% model error, i.e process noise
% it must be adjusted so we have a good SNR signal noise ratio
% s: signal noise
% w: process noise
W = eye(m)*(10*s);

% Observations are generated using the state space recursion.
% Function |dlmsmo| calculates the estimated states using Kalman smoother.
kfOut = dlmsmo(yTrain,F,V,x0,G,W,C0);
yEstm = kfOut.yhat;

% save smoothed estimate of x0
x0 = kfOut.x(:,1);


end