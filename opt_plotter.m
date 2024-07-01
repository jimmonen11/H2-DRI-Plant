close all

LMPcolor = '#4F7942';
CO2color = '#808080';
ndotcolor = "b";
metalcolor = "k";
IOcolor = "#A13D2D";

fontsize = 16;
fontsizetit = 21;
linewidthsize = 3;
hlinewidth = 3;

%time = out.metal.time/3600;


%%
% metal =  out.metal.data;
% m_IO = out.m_IO.data;
% LMP =  out.LMP.data;
% CO2 =  out.CO2.data;
% ndot =  out.ndot_H2O.data;

%%

%%
close all

load('params.mat')

t=0:1:8759;
% pandcdata = csvread('caiso_lmp_carbon_clean.csv', 1, 1);
% ndot_path = csvread('ndot1.csv');
% steel_path = csvread('steel1.csv');

pandcdata = csvread('ercot_lmp_carbon_clean.csv', 1, 1);
ndot_path = csvread('ndot_ERCOT.csv');


figure(1)
subplot(3,1,1)
box on
plot(t, pandcdata(:,1) ,  'linewidth', linewidthsize, 'color', LMPcolor )
ylabel('LMP ($/MWh)')
xticks(0:24:24*7)
xlim([0, 24*7])
yticks(0:25:100)
ylim([0, 125])


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

subplot(3,1,2)
box on
plot(t, pandcdata(:,2), 'linewidth', linewidthsize, 'color', CO2color )
ylabel('GSE (kg CO_2e/MWh)')
xticks(0:24:24*7)
xlim([0, 24*7])
%ylim([120, 320])
ylim([225, 445])
yticks(200:50:450)


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

cmap = colormap('cool');
x_cmap = linspace(0,1,length(cmap))';
alpha_vals = [0:0.125:1]';

for i =1:length(alpha_vals)
    colorCool(i,:) = interp1(x_cmap, cmap, alpha_vals(i));
end
subplot(3,1,3)
box on
hold on
for i = 1:9
  
    plot(t, ndot_path(:,i), 'linewidth', linewidthsize, 'color', colorCool(i,:) )
end

xlim([0, 24*7])
ylim([20, 170])
ylabel('Feed Water (mol/s)')
xticks(0:24:24*7)

H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)


figure(2)

cmap = colormap('cool');
x_cmap = linspace(0,1,length(cmap))';
alpha_vals = [0:0.125:1]';

for i =1:length(alpha_vals)
    colorCool(i,:) = interp1(x_cmap, cmap, alpha_vals(i));
end
subplot(2,1,2)
box on
hold on
for i = 1:9
    steel_path = ps(1)*ndot_path(:,i).^2 + ps(2)*ndot_path(:,i) + ps(3);
    plot(t, steel_path*1000, 'linewidth', linewidthsize, 'color', colorCool(i,:) )
end

xlim([0, 24*7])
ylim([5, 35])
ylabel('Liquid Steel (kg/s)')
xlabel('Time (hours)')
xticks(0:24:24*7)


pos = get(gca, 'position');
cbar_width = 0.02;  % Adjust this value as needed
cbar_height = pos(4);  % Match the height of the subplot
cbar_position = [pos(1) + pos(3) + 0.01, pos(2), cbar_width, cbar_height];
cbar = colorbar('location', 'SouthOutside', 'LineWidth',hlinewidth-1);
cbar.TickLength = 0.005; 
%cbar.Limits = [0:1];
cbar.YTick = 0:0.125:1;
cbar.FontSize = fontsize;
ylabel(cbar, '\alpha','FontWeight', 'bold','FontSize',fontsize+4)

pos1 = -0.3161;
pos2 = -16.9186;

pos3 = 168.0452;

text(pos1, pos2, 'Environmental', 'FontSize',fontsize, 'FontWeight', 'bold')
text(pos3, pos2, 'Economic', 'FontSize',fontsize, 'FontWeight', 'bold', 'HorizontalAlignment', 'right')


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

%%

figure(3)
subplot(3,1,1)
box on
plot(time, LMP ,  'linewidth', linewidthsize, 'color', LMPcolor )
ylabel('Elect. Price ($/MWh)')
xticks(0:24:24*7)
xlim([0, 24*7])
yticks(0:25:100)
ylim([-10, 85])


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

subplot(3,1,2)
box on
plot(time, CO2 , 'linewidth', linewidthsize, 'color', CO2color )
ylabel('Elect. Emissions (kg CO_2e/MWh)')
xticks(0:24:24*7)
xlim([0, 24*7])
yticks(100:50:300)


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

subplot(3,1,3)
box on
plot(time, ndotcolor, 'linewidth', linewidthsize, 'color', ndotcolor )
ylabel('Water Make-Up (mol/s)')
xlabel('Time (hrs)')
xticks(0:24:24*7)
xlim([0, 24*7])
ylim([0, 3.5])

H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)


%%

figure(3)

subplot(2,1,1)
box on
plot(time, metal*100, 'linewidth', linewidthsize, 'color', metalcolor )
ylabel('DRI Metallization (%)')
xticks(0:24:24*7)
xlim([0, 24*7])
ylim([93.5, 96.5])


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)


subplot(2,1,2)
box on
plot(time, m_IO , 'linewidth', linewidthsize, 'color', IOcolor )
ylabel('Iron Ore (kg/s)')
xlabel('Time (hrs)')
xticks(0:24:24*7)
xlim([0, 24*7])

H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)


%%

%%

figure(4)

subplot(3,1,3)
box on
plot(time, metal*100, 'linewidth', linewidthsize, 'color', metalcolor )
ylabel('DRI Metallization (%)')
xticks(0:24:24*7)
xlim([0, 24*7])
ylim([93.5, 96.5])
xlabel('Time (hours)')


H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)


subplot(3,1,1)
box on

plot(time, ndot*18/1000,  'linewidth', linewidthsize, 'color', ndotcolor )
ylabel('Feed Water (kg/s)')
ylim([0.5, 3])

xticks(0:24:24*7)
xlim([0, 24*7])

H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)

subplot(3,1,2)
plot(time, m_IO , 'linewidth', linewidthsize, 'color', IOcolor )
ylabel('Iron Ore (kg/s)')

xticks(0:24:24*7)
xlim([0, 24*7])

H = gca;
grid on
H.LineWidth = hlinewidth; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
%%

out.cost_total.data(end)/out.steel_total.data(end)
out.emissions_cost_total.data(end)/out.steel_total.data(end)/1000