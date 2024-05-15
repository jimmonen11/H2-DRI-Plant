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

time = out.metal.time/3600;


%%
metal =  out.metal.data;
m_IO = out.m_IO.data;
LMP =  out.LMP.data;
CO2 =  out.CO2.data;
ndot =  out.ndot_H2O.data;

%%

figure(1)
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

figure(2)

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

figure(3)

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