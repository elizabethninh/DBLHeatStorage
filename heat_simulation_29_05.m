clear all, close all, clc
%% Setup
% This section contains the parameters of the setup (start of cycle)

T_0 = 293;              %Starting temperature in [K], assumed for all fluids, gasses and materials
t = 0;                  %Time at the start [s]
t_final=1200;                 %Time at end cycle [s]
%% First order variables

%This section contains the properties and constants of the materials, but
%also fixed variables of the setup

%Copper tubing
r_inner_cu = 0.004;     %Inner radius copper tube [m]
r_outer_cu = 0.006;     %Outer radius copper tube [m]
roughness_cu = 0.0015;  %Pipe roughness [m]
length_cu = 6.85;       %Length copper tube [m]
epsilon_cu = 0.05;      %Emissivity copper [-]
rho_cu = 8960;          %Density copper [kg/m^3]
c_cu = 386;             %Specific heat copper [J/kg K]
k_cu = 400;             %Thermal conductivity copper [W/(m K)]

%PVC tubing
r_inner_pvc = 0.0528;   %Inner radius PVC tube [m]
r_outer_pvc = 0.055;    %Outer radius PVC tube [m]
length_pvc = 0.125;     %Length PVC tube [m]
epsilon_pvc = 0.9;      %Emissivity PVC [-]
rho_pvc = 1330;         %Density PVC [kg/m^3]
k_pvc = 0.19;           %Thermal conductivity PVC [W/(m K)]

%Polyurethane tubing
r_inner_pur = 0.004;    %Inner radius polyurethane tube [m]
r_outer_pur = 0.006;    %Outer radius polyurathane tube [m]
length_pur = 6;         %Length polyurethane tube [m]
epsilon_pur = 0.9;      %Emissivity polyurethane [-]
k_pur = 0.13;           %Thermal conducitivity polyurethane [W/(m K)]

%Aluminium plate
A_al = 1.085;           %Area aluminium plate [m^2]
epsilon_paint = 0.95;   %Emissivity black acrylic paint [-]
rho_al = 2700;          %Density aluminium [kg/m^3]
%Water
rho_water_20C = 998;        %Density water at 20C [kg/m^3]
c_water = 4148;         %Specific water [J/kg K]

%Air
rho_air = 1.29;         %Density air at 20C [kg/m^3]
c_air = 1006;           %Specific air [J/kg K]

%Other
sigma = 5.67*10^-8;     %Stefan Boltzmann constant [W/(m^2 K^4)]
flowrate = 3/60000;     %Flowrate pump [m^3/s]
E = 1000;               %Irradiance artificial sun [W/m^2]
%% Second order variables

%Copper tube
A_outer_cu = 2 * pi * r_outer_cu * length_cu;       %Outer surface area copper tube [m^2]
A_inner_cu = 2 * pi * r_inner_cu * length_cu;       %Inner surface area copper tube [m^2]
A_exposed_cu = A_outer_cu/2

V_cu = pi * r_inner_cu^2 * length_cu;               %Inner volume copper tube [m^3]
M_cu = V_cu * rho_cu

%Surface areas and volume PVC tube
A_outer_pvc = 2 * pi * r_outer_pvc * length_pvc;    %Outer area PVC tube [m^2]
A_inner_pvc = 2 * pi * r_inner_pvc * length_pvc;    %Inner area PVC tube [m^2]

V_pvc = pi * r_inner_pvc^2 * length_pvc;            %Inner volume PVC tube [m^3]

%Aluminium plate
V_al = A_al*0.002;                                  %Volume aluminium plate [m^3]
M_al = rho_al*V_al                                  %Mass aluminium plate [kg]

%Convective heat transfer coefficients
h_air = 10          %CVTH of still air [W/(m^2 K)], PLACEHOLDER
h_cu = 1700         %CVTH water in copper [W/(m^2 K)]
h_pur = 120.76      %CVTH water in polyurethane [W/(m^2 K)]
h_pvc = 4.21        %CVTH water in PVC [W/(m^2 K)]

%Other
V_system = V_cu+V_pvc;                  %Volume of system   

%% Plotting info

while t<t_final
   
    %Solar collector
    Q_rad_cu = E*length_cu*(r_outer_cu*2)*epsilon_paint;            %Heat addition radiation on copper tube [W]
    Q_rad_al = E*A_al*epsilon_paint;                                %Heat addition radiation on aluminium plate [W]
    Q_loss_conv_al = h_air*A_al*(T_0-T_final);                      %Heat loss convection aluminium plate [W]
    Q_loss_cond_al_cu = k_cu* A_exposed_cu*(T_al-T_cu)/(r_outer_cu*2-r_inner_cu*2) %dit klopt niet hlml maar weten A tussen plaat en buis niet %Heat loss conduction aluminium plate [W]
    Q_loss_cond_al =;       %Heat loss convection aluminium plate (other direction)
                                           
    R_sol_air = ;               %Conductive thermal resistance air 
    R_al =                      %Thermal resistance aluminium plate 
    R_sol_wood = ;              %Conductive thermal resistance wood setup
    
    Q_loss_rad_cu = ;       %Heat loss radiation copper tube [W]
    Q_loss_rad_al =;        %Heat loss radiation aluminium plate [W]
    
    T_al=T_al+(Q_rad_al-Q_loss_conv_al-Q_loss_cond_al-Q_loss_rad_al)/(m_al*c_al)  %Temperature of the aluminium plate [K]
    T_cu=T_cu+(Q_rad_cu-Q_loss_rad_cu-Q_loss_cond_al_cu)/(M_cu*c_cu)              %Temperature of the copper tube
    
    %Add heat flux copper tube
    
    %Polyurethane connection tube
    R_con_air = ;       %Convective thermal resistance air
    R_con_water = ;     %Convective thermal resistance water
    R_con_pur = ;       %Conductive thermal resistance polyurethane
    R_con_total = R_con_air + R_con_water + R_con_pur;     %Totral thermal resistance connection
   
    %Add heat flux polyurethane tubes
    
    %Storage vessel
    R_hsv_air =;        %Convective thermal resistance air
    R_hsv_water = ;     %Convective thermal resistance water
    R_hsv_pvc = ;       %Conductive thermal resistance PVC
    R_hsv_kingspan ;    %Conductive thermal resistance insulation Kingspan
    R_hsv_total = R_hsv_air + R_hsv_water + R_hsv_pvc + R_hsv_kingspan; 
    
    %Add heat flux Storage vessel
    
    %Heat transfer system
    Q_tot = ; %all fluxes added to each other   
    
    %Temperature water
    M_water = V_system*rho_water_20C;       %Volume of water inside system
    T_water = T_0+(Q_tot/M_water*c_water);      %Final temperature water [K]
    y(i) = T_water
    
    %Plotting steps code
    t = t + 1;          % 1 stands for 1 step for time
    x(i) = t;           % time is on x axis
    i = i + 1;          % 1 step added for plot
end
%% Actual plot

plot(x,y);
xlim([0, 1200]);
ylim([280, 350]);
xlabel('Time [s]');
ylabel('Temperature [K]');