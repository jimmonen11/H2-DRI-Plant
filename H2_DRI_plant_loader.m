%% Number of discretization modes and type of operation

n_furnace = 75;
H2only = true;

%% Furnace Geometry and Solids Parameters
% Geometry of furnace, pellets, etc. come from data from a plant in Quebec taken from 
% Hamadeh, H. (2017). Modélisation mathématique détaillée du procédé de réduction directe du minerai de fer. 1–143. https://theses.hal.science/tel-01740462v1/document
rho_p = 3528; % kg/m^3 density of a pellet
eps_bed = 0.5; % m^3 gas/m^3 total

r_furnace = 2.75; %m, radius of reducing section of furnace
h_furnace = 5; %m, height of reducing section of furnace
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
load('initcond_H2.mat')


x_CH4in = 0.0;
x_H2in = 0.90;
x_COin = 0.0;
x_H2Oin = 0.085;
x_CO2in = 0.0;
x_N2in = 0.015;
Gas_In_Flow = 2078.874*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s

T_gin = 925 + 273;

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
Kmetal = -0.0283;
tau_metal = 1.3*3600; 
tau_metal_desired = 60*3;


P_metal = tau_metal/(Kmetal*tau_metal_desired);
I_metal = P_metal/tau_metal;
D_metal = 0.5*I_metal;

Kd = 0.008729;
tau_d = 3600*(1.6);


tau_screw = 60*5; %5 minutes
K_screw = 0.5;

P_screw = 1/K_screw;
I_screw = P_screw/tau_screw;

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

CCsf = 250*ls_prod_kg*8760/1000; % Kruger from Bhaskar

CCelect = CC_elect*MWelect*1000; % $, cost of electrolyzers 

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

labor_cost = 19*DRI_prod_kg/3600/1000 + 53*ls_prod_kg/3600/1000; % $/s %Wood, Vogel cited by nature art
taxes = tax_rate*TIC/8760/3600;

CRF = 0.1019; %8% rate with 20 yr plant lifetime - from nature paper (r*(1+r)^20)/((1+r)^20-1)

%%

n=24;
day1 = 155;
day2 = 53;

%day1 = 1;
%day2 = 365;


pandcdata = csvread('caiso_lmp_carbon_clean.csv', 1, 1);
ndot_path = csvread('ndot1.csv');

ndot_path = ndot_path(day1*n+9:end,:);
pandcdata = pandcdata(day1*n+9:end,:);

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