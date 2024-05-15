function [mixed_profit, cost, emissions_cost, LCOS_bd, SCE_bd] = negative_profit( LCOE, elec_CO2rate, cIO, pLS, MWelect, MWeaf, MWcomp, MWpsa, MWcooltow, MWng, MWng_eaf,m_IO, m_steel, carbon_consume, lime_consume, ...
    CCeaf, CCelect, CCsf, CCfur, CCpsa, CCcomp, CCcooltow, CCbop, taxes, labor_cost, cNG, cCarbon, cLime, stack_replace,...
    cost_goal, emissions_goal)


% m_steel = tonne liquid steel per second

%==========================Variable O&M==============================


io_cost = m_IO*cIO/1000; % $/s, cost of iron ore
elec_cost = (MWelect+MWeaf+MWcomp+MWpsa+MWcooltow)*LCOE/3600; % $/s, cost of electricity
gas_cost = (MWng + MWng_eaf) *3.412*cNG/3600; % $/s, cost of ng

carbon_cost = carbon_consume*cCarbon; %$, cost of carbon
lime_cost = lime_consume*cLime; % cost of lime

stack_cost = MWelect*0.975*stack_replace/3600; % $/s %0.975 is rectifier efficiency - cost to replace stacks


%% LCOS
%CRF = 0.136;
CRF = 0.1019;

CCann = (CCeaf + CCelect + CCsf + CCfur + CCcomp + CCpsa + CCcooltow + CCbop)*CRF;
CC_cost = CCann/(8760*3600); % $/s for plant capital costs

LCOS = (io_cost + elec_cost + gas_cost + CC_cost + taxes + labor_cost + carbon_cost + lime_cost + stack_cost)/m_steel;
LCOS
cost = (io_cost + elec_cost + gas_cost + CC_cost + taxes + labor_cost + carbon_cost + lime_cost + stack_cost);
revenue = m_steel*cost_goal;

profit  = revenue - cost;

LCOS_bd = [CCelect*CRF/(8760*3600) CCsf*CRF/(8760*3600) CCeaf*CRF/(8760*3600) ...
    CCfur*CRF/(8760*3600)+CCcomp*CRF/(8760*3600)+CCpsa*CRF/(8760*3600)+CCcooltow*CRF/(8760*3600)  CCbop*CRF/(8760*3600)...
    elec_cost  gas_cost  io_cost   stack_cost  carbon_cost + lime_cost taxes/2 labor_cost taxes/2 ]/m_steel;

% ==========================CO2 Emissions==============================
% 
ng_CO2rate = 198; %kg CO2/MWh nat gas -from my 1st paper
ng_CO2rate_us = 47; %kg CO2e nat upstream - from Incentives for Clean Hydrogen Production in the Inflation Reduction Act

lime_CO2rate = 0.05; %kg CO2/kg lime consumed
carbon_CO2rate = 0.98; %kg CO2/kg carbon consumed
io_CO2rate = 0.02; %kg CO2/kg iron ore
pell_CO2rate = 0.16; %kg CO2/kg iron ore
og_CO2rate = 0.101; %kg CO2/kg ls
% 
elec_CO2 = (MWelect+MWeaf+MWcomp+MWpsa+MWcooltow)*elec_CO2rate/3600; %kg CO2/s from electricity
ng_CO2 = (MWng + MWng_eaf)*ng_CO2rate/3600; % kg CO2/s from natural gas combustion
 
EAFog_CO2 = og_CO2rate*m_steel*1000;  %kg CO2/s from combustion of EAF off gas
ng_CO2us = (MWng + MWng_eaf)*ng_CO2rate_us/3600; %kg CO2/s from upstream natural gas use


lime_CO2 = lime_CO2rate*lime_consume; %kg CO2/s from lime production
carbon_CO2 = carbon_CO2rate*carbon_consume; %kg CO2/s from coal/coke production
io_CO2 = io_CO2rate*m_IO; %kg CO2/yr from iron ore mining
pell_CO2 = pell_CO2rate*m_IO;  %kg CO2/yr from iron ore pelletizing
mineralCO2 = lime_CO2 + carbon_CO2 + io_CO2; %kg CO2/yr from mineral production/mining excluding pelletizing

SCE = (elec_CO2 + ng_CO2 + EAFog_CO2 + ng_CO2us + pell_CO2 + mineralCO2)/m_steel/1000 % t CO2e/ t ls
SCE_bd = [ng_CO2  EAFog_CO2  elec_CO2  ng_CO2us   pell_CO2  mineralCO2]/m_steel/1000; % t CO2e/ t ls

emissions_cost = (elec_CO2 + ng_CO2 + EAFog_CO2 + ng_CO2us + pell_CO2 + mineralCO2); %kg/s

%1.32

emissions_revenue = m_steel*emissions_goal*1000; %kg/s

emissions_profit = (emissions_revenue - emissions_cost);

%mixed_profit = -(profit + emissions_profit);


%LCOS
LCOSscaled = LCOS/cost_goal;
%SCE
%emissions_goal
SCEscaled = SCE/emissions_goal;

mixed_profit = (LCOSscaled + SCEscaled);
