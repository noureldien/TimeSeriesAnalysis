%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

return;

%Estimate an ARI model.
load iddata9 z9
y = z9.y;
model = ar(y, 4, 'ls');
compare(y,model,2);

return;

% %Use option set to choose 'ls' estimation approach
% % and to specify that covariance matrix should not be estimated.
% y = rand(100,1);
% opt = arOptions('Approach', 'ls', 'EstCovar', false);
% model = ar(y, N, opt);







