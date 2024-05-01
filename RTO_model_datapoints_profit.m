%Script to run profit and emissions profit model to create datapoints for RTO
%model. 

% clc; clear all; close all;
% 
n_samp = 1000; %number of samples to get points from the simulink model at
% n_dim = 5; %number of dimensions for the model, 2 decision vars
% n_match = 10000; %number of 'matches' for to create low discrepency grid
% % 
% winbox = minimax(n_samp, n_dim, n_match);
% 
% % 
% winbox(:,1) = scalefun(250,1000,winbox(:,1)); % molar flow rate of electrolyzer
% winbox(:,2) = scalefun(-15,250, winbox(:,2)); % LCOE
% winbox(:,3) = scalefun(50,360, winbox(:,3)); % kg/MWh from grid
% winbox(:,4) = scalefun(75,170, winbox(:,4)); % buying price of iron ore
% winbox(:,5) = scalefun(600,1700, winbox(:,5)); % sell price of liquid steel


%%

mixed_p = zeros(1,n_samp);

for i = 1:n_samp

    ndot_target = winbox(i,1);
    inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target);
    mixed_p(i) = negative_profit(winbox(i,2), winbox(i,3),winbox(i,4),winbox(i,5),...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
    CCeaf, CCelect, CCsf, CCbop, tax_rate, labor_cost, cNG, cCarbon, cLime, stack_replace);

end

p = polyfitn(winbox,mixed_p,2); %crete a quadratic model

polyvaln(p, [500, 25, 200, 130, 1000]) % test the model

%%
DV = 5; % number of decision vars
% 
% 
% diagonal = 2*[p.Coefficients(1) p.Coefficients(7) p.Coefficients(12) p.Coefficients(16) p.Coefficients(19)]
% 
% top = zeros(1, DV);
H = [  p.Coefficients(1)*2 p.Coefficients(2) p.Coefficients(3) p.Coefficients(4) p.Coefficients(5);...
       p.Coefficients(2)   p.Coefficients(7)*2 p.Coefficients(8) p.Coefficients(9) p.Coefficients(10);
       p.Coefficients(3)   p.Coefficients(8) p.Coefficients(12)*2 p.Coefficients(13) p.Coefficients(14);     
       p.Coefficients(4)   p.Coefficients(9) p.Coefficients(13) p.Coefficients(16)*2 p.Coefficients(17); 
       p.Coefficients(5)   p.Coefficients(10) p.Coefficients(14) p.Coefficients(17) p.Coefficients(19)*2;];

f = [p.Coefficients(6) p.Coefficients(11) p.Coefficients(15) p.Coefficients(18) p.Coefficients(20)]';



%%

lcoe = 75;
co2 = 400;
io = 125;
steel = 825;

% Aeq = diag([0, 1, 1 ,1, 1]);
% beq = [0, lcoe, co2 ,io, steel]';

Aeq = [];
beq = [];



% lb = [250, -15, 50, 75, 600]';
% ub = [1000, 250, 360, 170, 1700]';

lb = [250 lcoe co2 io steel];
ub = [1000 lcoe co2 io steel];

A = [];
b = [];

[x, p_mix] = quadprog(H,f,A,b,Aeq,beq,lb,ub)


%opt = -(p.Coefficients(2)*lcoe +  p.Coefficients(3)*co2 +  p.Coefficients(4)*io + p.Coefficients(5)*steel)/(2*p.Coefficients(1))
%%

lcoe = 59;
co2 = 329;
io = 160.96;
steel = 1100;

tp = 750;
for i = 1:tp
    in = [250+i lcoe co2 io steel];
    n_plot(i) = 250 + i;
    p_plot(i) = polyvaln(p, in);
end

close all
plot(n_plot, p_plot)