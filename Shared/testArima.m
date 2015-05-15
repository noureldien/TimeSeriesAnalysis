
model = arima(2,1,2);
[estModel, estParmCov, logL, info] = estimate(model, yTrain);
[E, V] = infer(estModel, yTrain);
%[Y,YMSE,V] = forecast(model, 48, 'X0', X, 'Y0', Y(100:length(Y)), 'XF', XF);

figure(1); clf;
hold on;
grid on;
box on;
plot(yTrain);
plot(yTrain+E);

return;