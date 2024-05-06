%Script to run profit and emissions profit model to create datapoints for RTO
%model. 

% clc; clear all; close all;
% 
n_samp = 1000; %number of samples to get points from the simulink model at
% n_dim = 3; %number of dimensions for the model, 2 decision vars
% n_match = 10000; %number of 'matches' for to create low discrepency grid
% % 
% 
% winbox = minimax(n_samp, n_dim, n_match);



%%

LCOEmin = -7;
LCOEmax = 2530;
CO2min = 135;
CO2max = 505;


winbox(:,1) = scalefun(31,144,winbox(:,1)); % molar flow rate of electrolyzer
winbox(:,2) = scalefun(LCOEmin,LCOEmax, winbox(:,2)); % LCOE
winbox(:,3) = scalefun(CO2min,CO2max, winbox(:,3)); % kg/MWh from grid

% winbox(:,4) = scalefun(75,170, winbox(:,4)); % buying price of iron ore
% winbox(:,5) = scalefun(600,1700, winbox(:,5)); % sell price of liquid steel


%%

io = 125;
steel = 825;

mixed_p = zeros(1,n_samp);

for i = 1:n_samp

    ndot_target = winbox(i,1);
    inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target);
    mixed_p(i) = negative_profit(winbox(i,2), winbox(i,3),io,steel,...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
    CCeaf, CCelect, CCsf, CCbop, tax_rate, labor_cost, cNG, cCarbon, cLime, stack_replace);

end

p = polyfitn(winbox,mixed_p,2) %crete a quadratic model

%polyvaln(p, [500, 25, 200, 130, 1000]) % test the model

%%
% DV = 5; % number of decision vars


% H = [  p.Coefficients(1)*2 p.Coefficients(2) p.Coefficients(3) p.Coefficients(4) p.Coefficients(5);...
%        p.Coefficients(2)   p.Coefficients(7)*2 p.Coefficients(8) p.Coefficients(9) p.Coefficients(10);
%        p.Coefficients(3)   p.Coefficients(8) p.Coefficients(12)*2 p.Coefficients(13) p.Coefficients(14);     
%        p.Coefficients(4)   p.Coefficients(9) p.Coefficients(13) p.Coefficients(16)*2 p.Coefficients(17); 
%        p.Coefficients(5)   p.Coefficients(10) p.Coefficients(14) p.Coefficients(17) p.Coefficients(19)*2;];
% 
% f = [p.Coefficients(6) p.Coefficients(11) p.Coefficients(15) p.Coefficients(18) p.Coefficients(20)]';


H = [  p.Coefficients(1)*2 p.Coefficients(2) p.Coefficients(3);
       p.Coefficients(2)   p.Coefficients(5)*2 p.Coefficients(6);
       p.Coefficients(3)   p.Coefficients(6) p.Coefficients(8)*2;];

f = [p.Coefficients(4) p.Coefficients(7) p.Coefficients(9)]';



%%
lcoe = 46;
co2 = 230;

% Aeq = diag([0, 1, 1 ,1, 1]);
% beq = [0, lcoe, co2 ,io, steel]';

Aeq = [];
beq = [];


% lb = [250, -15, 50, 75, 600]';
% ub = [1000, 250, 360, 170, 1700]';

lb = [35 lcoe co2 io steel];
ub = [140 lcoe co2 io steel];

A = [];
b = [];

[x, p_mix] = quadprog(H,f,A,b,Aeq,beq,lb,ub)


%opt = -(p.Coefficients(2)*lcoe +  p.Coefficients(3)*co2 +  p.Coefficients(4)*io + p.Coefficients(5)*steel)/(2*p.Coefficients(1))



tp = 139-35;
for i = 1:tp
    %in = [34+i lcoe co2 io steel];
    in = [34+i lcoe co2];
    n_plot(i) = 34 + i;
    p_plot(i) = polyvaln(p, in);
end

close all
plot(n_plot, p_plot)

% winbox(:,1) = scalefun(250,1000,winbox(:,1)); % molar flow rate of electrolyzer
% winbox(:,2) = scalefun(-15,250, winbox(:,2)); % LCOE
% winbox(:,3) = scalefun(50,360, winbox(:,3)); % kg/MWh from grid
% winbox(:,4) = scalefun(75,170, winbox(:,4)); % buying price of iron ore
% winbox(:,5) = scalefun(600,1700, winbox(:,5)); % sell price of liquid steel


ndot_target = x(1);
inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target+1);

mixed_p_scaler = negative_profit(lcoe, co2, io, steel,...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
    CCeaf, CCelect, CCsf, CCbop, tax_rate, labor_cost, cNG, cCarbon, cLime, stack_replace)

%%
options = optimoptions('quadprog','Display', 'off');

for i = 1:length(pandcdata)
    
    lb = [35 pandcdata(i,1) pandcdata(i,2)];
    ub = [140 pandcdata(i,1) pandcdata(i,2)];

    [x, p_mix] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[], options);
    
    ndot_opt(i) = x(1);

    if ndot_opt(i) > 140
        ndot_opt(i) = 35;

    end

    
    ndot_target = ndot_opt(i);
    inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target);

    mixed_p_opt(i) = negative_profit(pandcdata(i,1), pandcdata(i,2), io, steel,...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
    CCeaf, CCelect, CCsf, CCbop, tax_rate, labor_cost, cNG, cCarbon, cLime, stack_replace);


   
end

%% s

ind = 2029;

lb = [35 pandcdata(ind,1) pandcdata(ind,2)];
ub = [140 pandcdata(ind,1) pandcdata(ind,2)];

[x, p_mix] = quadprog(H,f,A,b,Aeq,beq,lb,ub)


%% s

