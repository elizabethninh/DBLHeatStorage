clear all, close all, clc
%% Setup
% This section contains the parameters of the setup (start of cycle)
i=1;                     %Value for plot
T_0 = 293;              %Starting temperature in [K], assumed for all fluids, gasses and materials
T_water = 293;              %Starting temperature water in [K]
T_al = 293;                 %Starting temperature aluminium in [K]
T_cu = 293;                 %Starting temperature copper in [K]
T_pvc = 293;                %Starting temperature PVC in [K]
T_pur = 293;                %Starting temperature polyurethane in [K]
T_air = 293;                %Starting temperature internal air in [K]
T_amb = 293;                %Starting temperature ambient (outside) air in [K]
t = 1;                      %Time at the start [s]
t_final=1200;               %Time at end cycle [s]

%% First order variables
%This section contains the properties and constants of the materials, but also fixed variables of the setup
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
hsv_pvc_d = 0.0022;     %Thickness of the HSV PCV [m]
r_outer_pvc = 0.055;    %Outer radius PVC tube [m]
r_inner_pvc = r_outer_pvc - hsv_pvc_d ;   %Inner radius PVC tube [m]
r_outer_double_pvc = r_outer_pvc + hsv_pvc_d; %Outer Radius of double-layered PVC tube (near endcaps) [m]
length_pvc = 0.130;     %Length/Height PVC tube [m] (this was shown to be slightly longer during first building phase with no real means to cut it down to 0.125m)
length_cap_pvc = 0.034; %Length/Height of a PVC Endcap [m]
epsilon_pvc = 0.9;      %Emissivity PVC [-]
rho_pvc = 1330;         %Density PVC [kg/m^3]
k_pvc = 0.19;           %Thermal conductivity PVC [W/(m K)]
d_air_endcap = 0.01;    %Distance between Kingspan Therma and the PVC endcap. Used to determine therma resistance
%Polyurethane tubing
r_inner_pur = 0.004;    %Inner radius polyurethane tube [m]
r_outer_pur = 0.006;    %Outer radius polyurathane tube [m]
d_pur = 0.002;          %Wall Thickness polyurthane tube [m]
length_pur = 6;         %Length polyurethane tube [m]
epsilon_pur = 0.9;      %Emissivity polyurethane [-]
k_pur = 0.13;           %Thermal conducitivity polyurethane [W/(m K)]
%Aluminium plate
A_al = 1.72*0.67;           %Area aluminium plate [m^2]
epsilon_paint = 0.95;   %Emissivity black acrylic paint [-] 
rho_al = 2700;          %Density aluminium [kg/m^3]
d_al = 0.002;           %thickness aluminium plate [m]
k_al = 237;             %Thermal conductivity aluminium [W/(m K)]
c_al = 900;             %Specific heat aluminium [J/kg K]
%Water
rho_water_20C = 998;    %Density water at 20C [kg/m^3]
c_water = 4148;         %Specific heat water [J/kg K]
%Air
rho_air = 1.29;         %Density air at 20C [kg/m^3]
c_air = 1006;           %Specific heat air [J/kg K]
k_air = 0.025;          %Thermal conductivity air [W/m K]
%General Setup Wood Properties
k_wood = 0.12;       %Thermal Conductivity of Pinewood [W/m k]                       
d_wood = 0.054;      %Wood border thcikness [m], simplification, in reality 2 mm air pocket within 
A_wood = 0.3107;     %Inner Surface Area wood layer [m^2], slight simplification
k_trespa = 0.03;     %Thermal Conductivity of Trespa bottom plate [W/m K]
d_trespa = 0.006;    %Thickness of trespa bottom [m]
%Kingspan Insulation Properties
ks_k = 0.022;     % Thermal Conductivity Kingspan insulation layer[W/(m*k)]
ks_d = 0.04;      % thickness of Kingspan insulation layer [m]
ks_spacer_d = 0.01; % Thickness of a single Kingspan spacer layer [m]
%Other
sigma = 5.67*10^-8;     %Stefan Boltzmann constant [W/(m^2 K^4)]
flowrate = 3/60000;     %Flowrate pump [m^3/s]
E = 1000;               %Irradiance artificial sun [W/m^2]
V_body_collector = 1.72*0.67*0.067; %Volume of the solar collector build area with glas plate [m^3]
%% Second order variables
%Copper tube
A_outer_cu = 2 * pi * r_outer_cu * length_cu;       %Outer surface area copper tube [m^2]
A_inner_cu = 2 * pi * r_inner_cu * length_cu;       %Inner surface area copper tube [m^2]
A_exposed_cu = A_outer_cu/2;
V_cu = pi * r_inner_cu^2 * length_cu;               %Inner volume copper tube [m^3]
M_cu = V_cu * rho_cu;
%Surface areas and volume PVC tube
A_outer_pvc = 2 * pi * r_outer_pvc * length_pvc;    %Outer area PVC tube [m^2]
A_inner_pvc = 2 * pi * r_inner_pvc * length_pvc;    %Inner area PVC tube [m^2]
V_pvc = pi * r_inner_pvc^2 * length_pvc;            %Inner volume PVC tube [m^3]
%Surface areas and volume Polyurethane Tubing
A_outer_pur = 2 * pi * r_outer_pur * length_pur;    %Outer area Polyurethane Tubing [m^2]
A_inner_pur = 2 * pi * r_inner_pur * length_pur;    %Inner area Polyurethane Tubing [m^2]
V_pur = pi * r_inner_pur^2 * length_pur;            %Inner volume Polyurethane Tubing [m^3]
%Aluminium plate
V_al = A_al*0.002;                                  %Volume aluminium plate [m^3]
M_al = rho_al*V_al;                                 %Mass aluminium plate [kg]
%The following values are lifted from the earlier matlab model and are used to calculate (average) surface area of the Kingspan Insulation.
hsv_al_width = 0.04017;                                                                % Width of Aluminum Tape reflector section, [m]
hsv_therma_width_in = (hsv_al_width + (2*(sind(18)*0.010)));                           %'width' of kingspan therma plate on the inside, [m]
hsv_therma_width_out = (hsv_al_width + (2*(sind(18)*0.010))) + 2*(tand(18)*ks_d);      %'width' of kingspan therma plate on the outside, [m]
hsv_therma_a_in = length_pvc * hsv_therma_width_in;                                    %Area of kingspan therma plate on inside, [m^2]
hsv_therma_a_out = length_pvc * (hsv_therma_width_out);                                %Area of kingspan therma plate on outside, [m^2]
hsv_therma_a_avg = (hsv_therma_a_out + hsv_therma_a_in)/2;                             %Surface area in middle of therma plate, [m^2]
%The Heat Storage Vessel is seperated into two parts: one for the parts with overlap of endcaps and one for the part without 
pvc_ec_inner_length = length_cap_pvc - hsv_pvc_d;                                                       %Inner length of PVC Endcap [m]
hsv_pvc_a_single = (length_cap_pvc - (pvc_ec_inner_length * 2)) * ((r_inner_pvc+r_outer_pvc)/2) * 2 * pi; %Average PVC Surface Area for Heat Storage Vessel section without overlap [m]
hsv_pvc_a_double = (pvc_ec_inner_length * 2) * ((pvc_ec_inner_length+r_outer_double_pvc)/2) * 2 * pi;         %Average PVC Surface Area for Heat Storage Vessel sections with overlap [m]
pvc_ec_a = (r_outer_pvc)^2 * pi;                                                                         %Surface area of a single endcap, [m^2]
% Air Insulation Values
hsv_air_avgdist1 = ((0.010+0.00795)/2)*((length_pvc - (pvc_ec_inner_length * 2))/length_pvc) + ((0.0078+0.00575)/2)*((pvc_ec_inner_length * 2)/length_pvc); %Average distance of the heat reflector to the outer wall of the PVC body [m], USES PRE CALCULATED VALUES.  Accounts for the double layer PVC near the endcap, [m]
hsv_air_avg_a1 = length_pvc * (hsv_al_width - (2*(sind(18)*0.005)));       %Surface area at the average point between pvc and reflector [m^2], (TODO implement avgdist 1)
hsv_air_avgdist2 = 0.00905;                                                %Average distance between reflector and kingspan therma plate [m]
hsv_air_avg_a2 = length_pvc * (hsv_al_width + (2*(sind(18)*0.005)));       %Surface are at the average point between reflector and kingspan therma [m^2], (TODO implement avgdist 2)
%Aluminum Reflector Tape Properties
d_al_t = 0.0001;                            % Thickness of an aluminum reflector segment [m]
a_al_t = hsv_al_width * length_pvc;         % Area of a single heat reflector element. 10 are present in total [m^2]
% Convective heat transfer coefficients
h_air = 5.0;         %CVTH of still air [W/(m^2 K)], Range (4.0 to 5.9)
h_water = 100;      %CVTH of water [W/(m^2 K)], Range (100 to 15000), Unused
h_pvc = 5;        %CVTH PVC [W/(m^2 K)], Range (5 to 9)
h_al = 7;        %CVTH aluminium [W/(m^2 K)], Range (7 to 10)
h_kingspan = 10;   %CTVH of Kingspan-Therma insulation [W/(m^2 K)], Range (10 to 30) (Aluminum surface)
h_pur = 5;        %CVTH of Polyurethane [W/(m^2 K)], PLACEHOLDER 
h_cu = 7;            %CVTH copper [W/(m^2 K)], Range (7 to 10) PLACEHOLDER
%Radiative heat transfer coefficients
h_r_pvc = 0.1;       %RDTH of pvc [W/(m^2 K)], PLACEHOLDER
h_r_al = 0.08;       %RDTH of aluminium-foil [W/(m^2 K)], PLACEHOLDER 
h_r_Kingspan = 0.08;  %RDTH of Kingspan-Therma insulation [W/(m^2 K)], PLACEHOLDER
h_r_pur = 0.5;       %RDTH of Polyurethane [W/(m^2 K)], PLACEHOLDER
%Overall Heat Transfer Coefficients, (may change due to change in CVTH)
U_cu = 1700;         %OHTC water in copper [W/(m^2 K)] TODO check this, has large influence on final temp
U_pur = 120.76;      %OHTC water in polyurethane [W/(m^2 K)]
U_pvc = 4.21;       %OHTC water in PVC [W/(m^2 K)]
%Other
V_system = V_cu+V_pvc+V_pur;                  %Volume of system  
A_contact_al_cu = length_cu * 0.002; %The contact patch area between the Aluminum and Copper. Simple lengthxwidth for area.
M_air = V_body_collector*rho_air;   
%Simplified without radiation
R_a = 1/((1/(h_pvc*(hsv_pvc_a_single+hsv_pvc_a_double))) + (log(r_outer_pvc/r_inner_pvc)/(2*pi*length_pvc*k_pvc))); %Convection from PVC into first air pocket, And conduction through PVC.
R_b = hsv_air_avgdist1 / (k_air * (10*hsv_air_avg_a1)); %Conduction through first air pocket
R_c = 1/((1/(h_al*(10*a_al_t))) + (d_al_t/k_al*(10*a_al_t))); %Convection from reflector to second air pocket, And conduction through reflector.
R_d = hsv_air_avgdist2 / (k_air * (10*hsv_air_avg_a2)); %Conduction through second air pocket
R_e = 1/((1/(h_kingspan*(10*hsv_therma_a_avg))) + (ks_d/(ks_k*(10*hsv_therma_a_avg)))); %Convection from Kingspan Therma to ambient, And conduction through Kingspan Therma.    
R_hsv_radial = R_a + R_b + R_c + R_d + R_e;
%Endcap (Assumed single endcap here. Second one will be accounted for in R_hsv_endcaps)
R_z = 1/((1/(k_pvc*pvc_ec_a)) + (hsv_pvc_d/(k_pvc*pvc_ec_a)));          %Convection from PVC into air gap, and conduction through pvc.
R_y = d_air_endcap/(k_air * pvc_ec_a);                                  %Conduction through air gap
R_x = 1/((1/(h_kingspan*pvc_ec_a)) + (ks_d/(ks_k*pvc_ec_a)));           %Convection from Kingspan Therma into outside atmosphere, and conduction through Kingspan Therma
R_hsv_endcaps = 1/(1/(R_z + R_y + R_x) + 1/(R_z + R_y + R_x));              %Total Equivalent Thermal resistance of both endcaps
R_hsv = 1/((1/R_hsv_radial) + (1/R_hsv_endcaps));
%misc thermal resistance 
R_al = d_al/(k_al * A_al);                    %Thermal resistance aluminium plate 
R_sol_air = 0.1;                                %Conductive thermal resistance air, PLACEHOLDER 
R_sol_wood = d_wood/(k_wood * A_wood) ;       %Conductive thermal resistance wood setup    
R_pur = 1/((1/(h_pur*(A_outer_pur))) + (log(r_outer_pur/r_inner_pur)/(2*pi*length_pur*k_pur))); %Thermal resistance of Polyurethane tube. Not accounting for rad
R_con = 1/((1/(h_cu*(0.2*(r_outer_cu*pi*2)))) + (log(r_outer_cu/r_inner_cu)/(2*pi*0.2*k_cu)))
R_trespa = d_trespa / (k_trespa * (A_al));
%Empty Arrays
T_air = repmat(0, 1,t_final);
T_cu = repmat(0, 1,t_final);
T_al = repmat(0, 1,t_final);
T_water = repmat(0, 1,t_final);
T_pur = repmat(0, 1,t_final);
T_0 = repmat(0, 1,t_final);
%Insertion of starting values (temporary)
T_0(1) = 293;                  %Starting temperature in [K], assumed for all fluids, gasses and materials
T_water(1) = 293;              %Starting temperature water in [K]
T_al(1) = 293;                 %Starting temperature aluminium in [K]
T_cu(1) = 293;                 %Starting temperature copper in [K]
T_pur(1) = 293;                %Starting temperature polyurethane in [K]
T_air(1) = 293;                %Starting temperature internal air in [K]
%% Plotting info
for i = 1:t_final
    %Water Mass
    %rho_water = (999.83953 + 16.945176 * (1.00024*T_water(i)) - 7.9870401*10^-3 * (1.00024*T_water(i))^3 - 46.17046*10^-6* (1.00024*T_water(i))^3 +105.56302*10^-9 * (1.00024*T_water(i))^4 - 280.54253*10^-12 * (1.00024*T_water(i))^5)/(1+16.897850*10^-3 * (1.00024*T_water(i)));   %Density water
    M_water = V_system*rho_water_20C;       %Volume of water inside system. Can produce weird values if the rest doesnt work, likely not cause of issues.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Addition of Heat energy through radiation, And loss of heat in collector solids and air - Needs improvement
    Qdot_rad_al = E*(A_al-A_exposed_cu)*epsilon_paint;                                %Heat addition radiation on aluminium plate [W]   
    Qdot_rad_cu = E*length_cu*(r_outer_cu* pi)*epsilon_paint;          %Heat addition radiation on copper tube [W] (only half of tube exposed directly to light)
    T_al(i) = T_al(i)+(Qdot_rad_al/(M_al*c_al));                       %Temperature of aluminum as it heats up from radiation [k]
    T_cu(i) = T_cu(i)+(Qdot_rad_cu/(M_cu*c_cu));                       %Temperature of aluminum as it heats up from radiation [k]    
    %Losses to inside of collector space
    Qdot_loss_conv_al(i) = h_air*(A_al-A_exposed_cu)*(T_al(i)-T_air(i));                      %Heat loss convection aluminium plate [W] 
    T_al(i) = T_al(i)-(Qdot_loss_conv_al(i)/(M_al*c_al));                       %Temperature of aluminum as it heats up from radiation [k]   
    Qdot_loss_conv_cu(i) = h_air*(A_outer_cu)*(T_cu(i)-T_air(i));
    T_cu(i) = T_cu(i)-(Qdot_loss_conv_cu(i)/(M_cu*c_cu)); 
    %Heating of air inside of collector
    T_air(i) = T_air(i) + ((Qdot_loss_conv_cu(i)+Qdot_loss_conv_cu(i))/(M_cu*c_cu)); 
    %Heated air going back into collector Al and Cu by means of convection (This stabilizes the internal temperatures and ensures that the ambient air's heat can eventually make its way into the liq)
    Qdot_gain_conv_al(i) = h_air*(A_al-A_exposed_cu)*(T_al(i)-T_air(i));     
    Qdot_gain_conv_cu(i) = h_air*(A_al-A_exposed_cu)*(T_al(i)-T_air(i)); 
    T_al(i) = T_al(i)+(Qdot_gain_conv_al(i)/(M_al*c_al));      
    T_cu(i) = T_cu(i)+(Qdot_gain_conv_cu(i)/(M_cu*c_cu));  
    T_air(i) = T_air(i) - ((Qdot_loss_conv_cu(i)+Qdot_loss_conv_cu(i))/(M_cu*c_cu)); 
    %Conduction from Aluminum to Copper
    Qdot_cond_al_to_cu(i) = h_cu*(A_contact_al_cu)*(T_al(i)-T_cu(i));
    T_cu(i) = T_cu(i)+(Qdot_cond_al_to_cu(i)/(M_cu*c_cu));  
    T_al(i) = T_al(i)-(Qdot_cond_al_to_cu(i)/(M_al*c_al)); 
    %----------------------------------------------------------------------
        %By this point no external losses in the collector are accounted for, following section transports away heat energy into the water.    
    Qdot_gain(i) = 1000*(A_inner_cu)*(T_cu(i)-T_water(i));                 %Fill sum of energy gains here [J/s or W]
    T_water(i) = T_water(i) + (Qdot_gain(i)/(M_water*c_water));            %New Temperature of water due to heat gain in a single second, right side is Temp change in said second
    T_cu(i) = T_cu(i)-(Qdot_gain(i)/(M_cu*c_cu));       
        %Aluminum Plate to bottom (assumes perfect contact)
    Qdot_loss_trespa(i) = (T_al(i)-T_amb)/R_trespa;                        %Heat loss through the trespa plate at the bottom
    T_al(i) = T_al(i)-(Qdot_loss_trespa(i)/(M_al*c_al)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Liquid - loss - This works as intended
    Qdot_loss_hsv(i) = (T_water(i)-T_amb)/R_hsv ;                           %Heat loss hsv
    Qdot_loss_pur(i) = (T_water(i)-T_amb)/R_pur;                            %Heat loss polyurethane tubing 
    Qdot_loss_con(i) = (T_water(i)-T_amb)/R_con;                            %Heat loss of the two 10cm long press-fit coupling tube sections sticking out of HSV
    Qdot_loss(i) = Qdot_loss_hsv(i) + Qdot_loss_pur(i)+Qdot_loss_con(i);    %Sum all of the heat losses per second instance here
    T_water(i) = T_water(i)-(Qdot_loss(i)/(M_water*c_water));               %New Temperature of water due to heat loss in a single second, right side is Temp change in said second
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Bringing the values that were determined in this time step to the next one
    T_water(i+1) = T_water(i);                                              %brings current water temp over to the next time increment
    T_al(i+1) = T_al(i);
    T_cu(i+1) = T_cu(i);
    T_air(i+1) = T_air(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plotting code
    y(i) = T_water(i);                                                      %Inserts current temp into its respective place in y for plotting
    t = t + 1;                                                              %1 stands for 1 step for time
    x(i) = i;                                                               %time is on x axis
    i = i + 1;                                                              %1 step added for plot
end
%% Actual plot
%Temperature
figure(1)
plot(x,y)
hold on
plot(x(1:60:t_final),y(1:60:t_final), '.', 'Markersize', 5)
hold off
xlim([0, t_final]);
ylim([293, 343]);
xlabel('Time [s]');
ylabel('Temperature [K]');
%Plot of Heat Energy Gain per Second
figure(2)
plot(x,Qdot_gain)
xlim([0, t_final]);
ylim([0, 800]);
xlabel('Time [s]');
ylabel('Heat Energy Gain in a second [J/s & W]');
%Plot of Heat Energy Loss per Second
figure(3)
plot(x,Qdot_loss)
xlim([0, t_final]);
ylim([0, 800]);
xlabel('Time [s]');
ylabel('Heat Energy Loss in a second [J/s & W]');

    %Q_rad_cu = E*length_cu*(r_outer_cu*2 * pi)*epsilon_paint;            %Heat addition radiation on copper tube [W]
    %Q_rad_al = E*A_al*epsilon_paint;                                %Heat addition radiation on aluminium plate [W]
    %Q_loss_conv_al(i) = h_air*A_al*(T_al(i)-T_air(i));                      %Heat loss convection aluminium plate [W]
    %Q_loss_cond_al_cu = k_cu* A_exposed_cu*(T_al(i)-T_cu(i))/(r_outer_cu*2-r_inner_cu*2);  %dit klopt niet hlml maar weten A tussen plaat en buis niet %Heat loss conduction aluminium plate [W]
    %Q_losscond(i)= Q_loss_cond_al_cu;
    
    %Q_loss_cond_al(i) = (T_al(i)-T_air(i))/R_al ;                                       %Heat loss convection aluminium plate (other direction)
    %Q_loss_rad_cu(i) = sigma * epsilon_paint * A_outer_cu* (T_cu(i)^4 - T_air(i)^4);    %Heat loss radiation copper tube [W]
    %Q_loss_rad_al(i) = sigma * epsilon_paint * A_al * (T_al(i)^4- T_air(i)^4);          %Heat loss radiation aluminium plate [W]