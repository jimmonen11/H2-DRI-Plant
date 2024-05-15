%% Number of discretization modes and type of operation

n_furnace = 75;

changenodes = false; %set to true if plan on changing n_furnace to something other than 75

H2only = true; % set to true for a furnace that only uses H2

steady_state = true; % set to true if want constant inlet conditions to be constant

% set one of these to true to run a dynamic case - only can use 1
% stepcase 1 has to have H2 false
stepcase1 = false;
stepcase2 = false;
stepcase3 = false;

if steady_state+stepcase1+stepcase2+stepcase3 > 1
    error('Trying to run a steady state case and a step change case.')
end

if H2only + stepcase1 > 1
    error('Hydrogen only case cannot be used with step case 1.')
end

%% Furnace Geometry and Solids Parameters
% Geometry of furnace, pellets, etc. come from data from a plant in Quebec taken from 
% Hamadeh, H. (2017). Modélisation mathématique détaillée du procédé de réduction directe du minerai de fer. 1–143. https://theses.hal.science/tel-01740462v1/document
rho_p = 3528; % kg/m^3 density of a pellet
eps_bed = 0.5; % m^3 gas/m^3 total

r_furnace = 2.75; %m, radius of reducing section of furnace
h_furnace = 10; %m, height of reducing section of furnace
r_p = 15e-3/2; %m, radius of iron ore pellets

%% Constants

% Molar mass of species in kg/mol
MM_Fe2O3 = 159.69/1000;
MM_Fe3O4 = 231.53/1000;
MM_FeO = 71.844/1000;
MM_Fe = 55.845/1000;
MM_C = 12.011/1000; 
MM_Gan = 60.08/1000;

% Molar mass of species in g/mol
MM_H2O = 18.015;
MM_H2 = 2.016;
MM_N2 = 28.013;
MM_CO = 28.01;
MM_CO2 = 44.01;
MM_CH4 = 16.04; 

% gas diffusion properties - MM and diffusion volumes
H2O_diff_prop = [MM_H2O, 13.1];
H2_diff_prop = [MM_H2, 2.31*2];

CO_diff_prop = [MM_CO, 18.9];
CO2_diff_prop = [MM_CO2, 26.9];

% ideal gas constant
R = 8.314; %(m^3*Pa)/(K*mol)

%% Solid Inlet Conditions
Solids_In_Flow = 44.1; %kg/s, flow of solids in

T_sin = -10 + 273; %K, solids temperature in
T_sin = 10 + 273; %K, solids temperature in

P_gin = 101325*1.8; %Pa, pressure of gas in (found via trial and error)

% weight fraction of pellets in
w_Fe2O3in = 0.9665;
w_Fe3O4in = 0;
w_FeOin = 0;
w_Fein = 0;
w_Cin = 0;
w_Ganin = 0.0335;

%% Gas Inlet Conditions

if H2only == true
    load('initcond_H2.mat')
    %load('initcond_H2_10m.mat') % can load 10 m results if want but change h_furnace
    h_furnace = 5;
    
    x_CH4in = 0.0;
    x_H2in = 0.962059;
    x_COin = 0.0;
    x_H2Oin = 0.019314;
    x_CO2in = 0.0;
    x_N2in = 0.018627;
    Gas_In_Flow = 2078.874*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s
    T_gin = 947 + 273;
    
    ndotin = 2078.874;
    
else
    load('initcond_NG.mat')
    
    % Steady state NG-DRI conditions
    x_CH4in = 0.105824828;
    x_H2in = 0.489526511;
    x_COin = 0.320711862;
    x_H2Oin = 0.042557156;
    x_CO2in = 0.024;
    x_N2in = 0.017379642;
    Gas_In_Flow = 2228.1*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000; %kg/s
    T_gin = 947 + 273;


end

if changenodes == true %this 'helps' change nodes, but will still required sim time to sort itself out

    T_ginit = interp1([1:1:length(T_ginit)],T_ginit, linspace(1,length(T_ginit),n_furnace));
    T_sinit = interp1([1:1:length(T_sinit)],T_sinit, linspace(1,length(T_sinit),n_furnace));
    
    c_H2Oinit = interp1([1:1:length(c_H2Oinit)],c_H2Oinit, linspace(1,length(c_H2Oinit),n_furnace));
    c_H2init = interp1([1:1:length(c_H2init)],c_H2init, linspace(1,length(c_H2init),n_furnace));
    c_N2init = interp1([1:1:length(c_N2init)],c_N2init, linspace(1,length(c_N2init),n_furnace));
    c_COinit = interp1([1:1:length(c_COinit)],c_COinit, linspace(1,length(c_COinit),n_furnace));
    c_CO2init =interp1([1:1:length(c_CO2init)],c_CO2init, linspace(1,length(c_CO2init),n_furnace));
    c_CH4init =interp1([1:1:length(c_CH4init)],c_CH4init, linspace(1,length(c_CH4init),n_furnace));
    
    c_Feinit = interp1([1:1:length(c_Feinit)],c_Feinit, linspace(1,length(c_Feinit),n_furnace));
    c_FeOinit = interp1([1:1:length(c_FeOinit)],c_FeOinit, linspace(1,length(c_FeOinit),n_furnace));
    c_Fe3O4init = interp1([1:1:length(c_Fe3O4init)],c_Fe3O4init, linspace(1,length(c_Fe3O4init),n_furnace));
    c_Fe2O3init = interp1([1:1:length(c_Fe2O3init)],c_Fe2O3init, linspace(1,length(c_Fe2O3init),n_furnace));
    c_Cinit = interp1([1:1:length(c_Cinit)],c_Cinit, linspace(1,length(c_Cinit),n_furnace));
    
    
    nr1init = interp1([1:1:length(nr1init)],nr1init, linspace(1,length(nr1init),n_furnace));
    nr2init = interp1([1:1:length(nr2init)],nr2init, linspace(1,length(nr2init),n_furnace));
    nr3init = interp1([1:1:length(nr3init)],nr3init, linspace(1,length(nr3init),n_furnace));
    ndotinit = interp1([1:1:length(ndotinit)],ndotinit, linspace(1,length(ndotinit),n_furnace));

end

x_sumin = x_CH4in  + x_H2in + x_COin +x_H2Oin + x_CO2in + x_N2in; %check to make sure equal to 1

if steady_state == true
    
    %inlet values are the same throughout simulation
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*1;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin;

end

if stepcase1 == true

    % 25% turndown of syngas with CH4/N2
    x_CH4step = 0.11138;
    x_H2step = 0.76842;
    x_COstep = 0.08439;
    x_H2Ostep = 0.01120;
    x_CO2step = 0.00632;
    x_N2step = 0.01829;
    Gas_In_Flow_Step = 14.16536754;
    T_ginstep = 947+273;

end


if stepcase2 == true 
    
    %Move to 80% of flow
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*0.8;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin;

end

if stepcase3 == true 
    
    %Move to 80% of flow
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*1;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin-100;

end

x_sumstep = x_CH4step  + x_H2step + x_COstep +x_H2Ostep + x_CO2step + x_N2step; %check to make sure equal to 1

x_CH4in = 0.0;
x_H2in = 0.90;
x_COin = 0.0;
x_H2Oin = 0.085;
x_CO2in = 0.0;
x_N2in = 0.015;
Gas_In_Flow = 2078.874*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s

T_gin = 950 + 273;

P_g_sp = 101325*2.5;

ndotin = 2706;

%% Furnace Geometry Again - Need to Know Height of Furnace

A_furnace = pi*r_furnace^2; %m^2, c.s. area of flow
A_furnace_pel = A_furnace*(1-eps_bed); %c.s. area for pellet flow - excludes gas lanes
A_furnace_gas = A_furnace*(eps_bed); %c.s. area for gas flow - excludes solids lanes

V_p = 4/3*pi*r_p^3; %m^3, volume of a pellet
V_furnace = A_furnace*h_furnace; %m^3, volume of reducing section of furnace
V_pellet_bed = ((4/3)*pi*r_p^3)/(1-eps_bed); %m^3, volume of pellets in furnace

dz = (h_furnace/(n_furnace-1)); %m, spacing of nodes

n_pellets = V_furnace/V_pellet_bed; % no. of pellets in reducing section
n_pellets_dz = n_pellets*dz/h_furnace; % no. of pellets per node

% Volume of gas and solid in each spatial node - constant
V_g = (4/3)*pi*r_p^3*n_pellets_dz*(eps_bed/(1-eps_bed)); % m^3, volume of gas
V_s = (4/3)*pi*r_p^3*n_pellets_dz; % m^3, volume of solid

a_b = 6*(1-eps_bed)/(r_p*2); %m^2/m^3, suface area for gas solid heat exchange, Wagner thesis
A_wall = (2*pi*r_furnace*dz); %m^2 surface area of furnace

tau = 1; % seconds, time constant for 1st order step changes

%% N2 Calculation

Vdot_s = Solids_In_Flow/rho_p;
Vdot_tot = Vdot_s/(1-eps_bed);
Vdot_g = Vdot_tot*eps_bed;

P_amb = 101325; 

ndot_inert = (P_amb*Vdot_g)/(T_sin*R);
%% SOE electrolyzer

SOEeff = 37.5; %kWh/kg H2

%% PID tuning

%from step test
Kmetal = -0.0435;
tau_metal = 3.2*3600; 

% Kd = 5.29e-4;
% tau_d = 3600*2;

%Kd = 3.2317e-4;
% Kd = 38.96;

%Kd = 0.1162;
Kd = 0.00728;
tau_d = 5200;
tau_d = 3600;


tau_metal_desired = 3600*0.15;
tau_metal_desired = 3600*0.1;

P_metal = tau_metal/(Kmetal*tau_metal_desired)
I_metal = P_metal/tau_metal
D_metal = 0.5*I_metal;

% Ktemp = 3.2e-4;
% tau_temp = 3600;
% 
% tau_temp_desired = 3600*0.5;
% 
% P_temp = tau_temp/(Kmetal*tau_temp_desired);
% I_temp = P_metal/tau_temp;

%% Reycle Compressor 

% Polytropic coefficient assumed to be nearly constant at 1.4

eta_isen = 0.75; % isentropic efficiency of compressor
eta_motor = 0.95; % compressor motor efficiency

%% Purge and Combustion

eta_comb = 0.8; %combustion efficiency


%% EAF

EAF_elec_req = 375; % kWhe/tonne steel
m_O2rate = 45.072; % kg O2/tonne steel
oxy_eff_tonne = 0.803; % MWhe/tonne O2


lime_CO2rate = 0.05; % kg CO2/kg lime consumed
carbon_CO2rate = 0.98; % kg CO2/kg carbon consumed
io_CO2rate = 0.02; % kg CO2/kg iron ore
pell_CO2rate = 0.16; % kg CO2/kg iron ore
og_CO2rate = 0.105; % kg CO2/kg ls -0.126 from supp of good journal, 0.105 from me

DRIsteel_conv = 1.075; %ratio of tonnes DRI to steel produced

m_steel = 30.07/1000; %tonne/s steel
ls_prod_kg = 30.07*3600; % kg/hr liquid steel
DRI_prod_kg = 32.32*3600; % kg/hr DRI
IO_input_kg = 45.07*3600; % kg/hr IO

m_O2 = m_O2rate*m_steel*3600; %kg/h O2 req
Qfurnace = 40.31; % MW
Pcomp = 10.77; % MW
MWelect = 240.452; %MW of electrolyzer required
Qwater = 6519.8; % m^3/hr for eqivalent 


%% Variable O&M

cIO = 130/1000; % $/kg, cost of iron ore 
cNG = 7; % $/MMBtu
cCarbon = 180/1000; % $/kg carbon
cLime = 100/1000; % $/kg lime
stack_replace = 3.25; % $/MWh-DC
MWhNGeaf = 0.08792; %MWh natural gas/tonne ls for EAF

hot_standby_eff = 0.009; %

%% Captial Costs

CC_elect = 2000; % $/kW - estimate that down from 2500

%CCsf = 49080*(IO_input_kg)^(0.6538)
CCsf = 250*ls_prod_kg*8760/1000; % Kruger from Bhaskar

CCelect = CC_elect*MWelect*1000; % $, cost of electrolyzers 

%CCeaf = 1132370*ls_prod_kg^0.4560 %$, cost of electric arc furnace
CCeaf = 160*ls_prod_kg*8760/1000; % Vogl from Bhaskar

CCpsa = 30622*m_O2^0.6357; % from Rosner all below
CCfur = 228860*Qfurnace^0.7848;
CCcomp = 6151202*Pcomp^0.71;
CCcooltow = 2*60812*Qwater^0.6303; %twice to account for harder methods...
CCbop = 69819*ls_prod_kg^0.5584 + 6320*ls_prod_kg^.8000 + 174548*ls_prod_kg^0.5583;

TIC = CCsf + CCelect + CCeaf + CCpsa+ CCfur+ CCcomp + CCcooltow + CCbop;

%% Taxes/Insurance and Labor

tax_rate = 0.02 + 0.02; % tax rate of annual capital cost also insurance and maintenance!
    % Wood gives another 2% for maintenance

%labor = 60*2080*62 %/yr 60 people making $62/hr

% labor = (51*40.85 + 93*30)*2080;
% 
% GA = 0.2*labor; % $/yr, 20% of labor
% 
% labor_cost = (labor + GA)/8760/3600 % $/s for labor

labor_cost = 19*DRI_prod_kg/3600/1000 + 53*ls_prod_kg/3600/1000; % $/s %Wood, Vogel cited by nature art


taxes = tax_rate*TIC/8760/3600;

CRF = 0.1019; %8% rate with 20 yr plant lifetime - from nature paper (r*(1+r)^20)/((1+r)^20-1)

% check this!
%BOPother = 69819*ls_prod_kg^0.5584 + 6320*ls_prod_kg^.8000 + 174548*ls_prod_kg^0.5583;
%cool_tower = 60812*(1907.7)^0.6303;

%CCbop = BOPother + cool_tower;

% 
% cooling_tower = 60812*(1907.7*CF)**0.6303
% BOPboil = MWboil*513*1000 #513 from my work
% 
% BOPother = 69819*ls_prod_kg**0.5584 + 6320*ls_prod_kg**.8000 + 174548*ls_prod_kg**0.5583

%%

% LMP = 57.13;
% GSE = 250.72;
% 
% [mixed_p, cost, emissions_cost, LCOS_bd, SCE_bd] = negative_profit(LMP, GSE, 130, 900,...
%     out.obj_fun_inputs.data(end,2), out.obj_fun_inputs.data(end,3), out.obj_fun_inputs.data(end,4), out.obj_fun_inputs.data(end,5), out.obj_fun_inputs.data(end,6), out.obj_fun_inputs.data(end,7), out.obj_fun_inputs.data(end,8), out.obj_fun_inputs.data(end,9), out.obj_fun_inputs.data(end,10), out.obj_fun_inputs.data(end,11), out.obj_fun_inputs.data(end,12),...
%     CCeaf, CCelect, CCsf, CCfur, CCpsa, CCcomp, CCcooltow, CCbop, taxes, labor_cost, cNG, cCarbon, cLime, stack_replace,...
%     10, 10)


%%

n=24;
day1 = 46;
day2 = 53;

%day1 = 1;
%day2 = 365;


pandcdata = csvread('caiso_lmp_carbon_clean.csv', 1, 1);
ndot_path = csvread('ndot1.csv');

% ndot_path = ndot_path(day1*n:day2*n,2);
% pandcdata = pandcdata(day1*n:day2*n,:);

ndot_path = ndot_path(:,3);

tlen = length(ndot_path);
t = 0:3600:3600*(tlen-1);
t = t';

runner = false;

if runner == true
    options = simset('SrcWorkspace','current'); %set the workspace to the current one
    output = sim('SOE_H2_DRI_plant',3600*24*6, options);
    
    obj_fun_inputs = output.obj_fun_inputs.data(end,:);
end