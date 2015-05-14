%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

%confband(t,out.x(1,:)',out.xstd(:,1));
%errorbar(t,y,2*V,'ok-');

% read csv and save it
% stockIndex = csvread('data/stockIndex.csv');
% stockIndex = flipud(stockIndex);
% save('data/stockIndex.mat', 'stockIndex');
% return;
% stockIndexMonthly = csvread('data/stockIndexMonthly.csv');
% stockIndexMonthly = flipud(stockIndexMonthly);
% save('data/stockIndexMonthly.mat', 'stockIndexMonthly');
% return;

% % read price to earning ratio (PER), flip it and save it as .mat
% dataBER = csvread('data/dataBER.csv');
% dataBER = flipud(dataBER);
% save('data/dataBER.mat', 'dataBER');
% return;

% % read the oil price and save it as .mat
% dataOIL = csvread('data/dataOIL.csv');
% save('data/dataOIL.mat', 'dataOIL');
% return;

% % read the NAPM and save it as .mat
% % ISM Manufacturing: PMI Composite Index (NAPM)
% dataPMI = csvread('data/dataPMI.csv');
% save('data/dataPMI.mat', 'dataPMI');
% return;

% % read the income and save it as .mat
% % ISM Manufacturing: PMI Composite Index (NAPM)
% dataINCOME = csvread('data/dataINCOME.csv');
% save('data/dataINCOME.mat', 'dataINCOME');
% return;

% % read the population and save it as .mat
% dataPOP = csvread('data/dataPOP.csv');
% save('data/dataPOP.mat', 'dataPOP');
% return;

% % read the unemployment and save it as .mat
% dataUNRATE = csvread('data/dataUNRATE.csv');
% save('data/dataUNRATE.mat', 'dataUNRATE');
% return;

% corporate profit is a bit annoying, because
% it need to be changed expanded from quarterly to monthly
% dataPROFIT1 = csvread('data/dataPROFIT.csv');
% dataPROFIT2 = zeros(size(dataPROFIT1, 1) * 3, 1);
% for i=1:length(dataPROFIT1)
%     dataPROFIT2( 1 + ((i-1)*3) : (i*3) ) = dataPROFIT1(i);
% end
% dataPROFIT = dataPROFIT2(5:end-1);
% save('data/dataPROFIT.mat', 'dataPROFIT');

% % now, form the big matrix that comprise all the data
% load('data/dataBER');
% load('data/dataOIL');
% load('data/dataPMI');
% load('data/dataINCOME');
% load('data/dataPROFIT');
% load('data/dataPOP');
% load('data/dataUNRATE');
% dataFeatures = [dataBER, dataOIL, dataPMI, dataINCOME, ...
%                 dataPROFIT, dataPOP, dataUNRATE];
% save('data/dataFeatures.mat', 'dataFeatures');












