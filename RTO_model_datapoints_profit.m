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

load('winbox_noscale.mat')

LCOEmin = -10;
LCOEmax = 143;
CO2min = 60;
CO2max = 390;

winbox(:,1) = scalefun(35.1,149.9,winbox(:,1)); % molar flow rate of electrolyzer
winbox(:,2) = scalefun(LCOEmin,LCOEmax, winbox(:,2)); % LCOE
winbox(:,3) = scalefun(CO2min,CO2max, winbox(:,3)); % kg/MWh from grid

%%

io = 130;
steel = 825;

cost_goal = 705.3;
emissions_goal = 1.103;

% cost_goal = 650;
% emissions_goal = 1;

mixed_p = zeros(1,n_samp);

for i = 1:n_samp

    ndot_target = winbox(i,1);
    inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target);
    
    [mixed_p(i), cost(i), emissions_cost(i), cost_noelect(i), em_noelect(i), ~, ~] = negative_profit(winbox(i,2), winbox(i,3),io,steel,...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
    CCeaf, CCelect, CCsf, CCfur, CCpsa, CCcomp, CCcooltow, CCbop, taxes, labor_cost, cNG, cCarbon, cLime, stack_replace,...
    cost_goal, emissions_goal);

end

cost_x = [winbox(:,1) winbox(:,2)];
p_cost = polyfitn(cost_x,cost,2) %crete a quadratic model
polyvaln(p_cost, [150, 83])

cost_x = [winbox(:,1)];
p_cost_noelect = polyfitn(cost_x, cost_noelect,2) %crete a quadratic model
polyvaln(p_cost_noelect, [150])

emissions_cost_x = [winbox(:,1) winbox(:,3)];
p_emissions_cost = polyfitn(emissions_cost_x,emissions_cost,2) %crete a quadratic model
polyvaln(p_emissions_cost, [150, 83])

steel_x = obj_function_inputs(:,1); % all the ndot values
p_steel = polyfitn(steel_x,obj_function_inputs(:,10),2)

MW_x = obj_function_inputs(:,1); % all the ndot values
p_MW = polyfitn(MW_x, sum(obj_function_inputs(:,2:5),2),2)


cost_x = [winbox(:,1)];
p_em_noelect = polyfitn(cost_x, em_noelect,2)

%%

% H = [  p_cost.Coefficients(1)*2 p_cost.Coefficients(2) p_cost.Coefficients(3);
%        p_cost.Coefficients(2)   p_cost.Coefficients(5)*2 p_cost.Coefficients(6);
%        p_cost.Coefficients(3)   p_cost.Coefficients(6) p_cost.Coefficients(8)*2;];
% 
% f = [p_cost.Coefficients(4) p_cost.Coefficients(7) p_cost.Coefficients(9)]';


%%
lcoe = 100;
co2 = 400;

% Aeq = diag([0, 1, 1 ,1, 1]);
% beq = [0, lcoe, co2 ,io, steel]';

Aeq = [];
beq = [];


% lb = [250, -15, 50, 75, 600]';
% ub = [1000, 250, 360, 170, 1700]';

lb = [37.5 lcoe co2 io steel];
ub = [150 lcoe co2 io steel];

A = [];
b = [];

[x, p_mix] = quadprog(H,f,A,b,Aeq,beq,lb,ub)


%opt = -(p.Coefficients(2)*lcoe +  p.Coefficients(3)*co2 +  p.Coefficients(4)*io + p.Coefficients(5)*steel)/(2*p.Coefficients(1))



tp = 150-37;
for i = 1:tp
    %in = [34+i lcoe co2 io steel];
    in = [36 + i lcoe co2];
    n_plot(i) = 36 + i;
    p_plot(i) = polyvaln(p, in);
end

close all
plot(n_plot, p_plot)


%%
ndot_target = x(1);
inputs_interp = interp1q(obj_function_inputs(:,1), obj_function_inputs(:,2:end), ndot_target+1);

[mixed_p_scaler, LCOS, SCE] = negative_profit(lcoe, co2, io, steel,...
    inputs_interp(1), inputs_interp(2), inputs_interp(3), inputs_interp(4), inputs_interp(5), inputs_interp(6), inputs_interp(7), inputs_interp(8), inputs_interp(9), inputs_interp(10), inputs_interp(11),...
     CCeaf, CCelect, CCsf, CCfur, CCpsa, CCcomp, CCbop, taxes, labor_cost, cNG, cCarbon, cLime, stack_replace,...
     cost_goal, emissions_goal)

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

