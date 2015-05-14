
% generate data
t = (1:1:(6*12))';
n = length(t);
% obs error std
s = 0.1;
y = 2 + 0.5*(t-0)/12.*((t-0)<2*12) + 0.5*(2*12)/12.*((t-0)>=2*12)  + 0.8*sin((t-0)/12*2*pi) + randn(n,1)*s;

% fit model
w123 = [0.0000,3.1e-5,1e-8,1e-8,1e-8,1e-8];
clear options
options.ns = 12;
options.trig = 2;
options.opt = 1;
out = dlmfit(y,s,w123,[],[],[],options);

% plot result
figure(1); clf
hold on;
grid on;
box on;
plot(t,y,'b');
plot(t,out.yhat,'r');
xlabel('time');
ylabel('y');
title('observations');










