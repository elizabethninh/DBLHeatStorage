clear all, close all, clc
%% This Version only takes the General Setup (PolyUrethane Tubes) and the Solar Collector (SSA8 version) into account. Setup values are coppied from earlier model

%% Setup
% This section contains the parameters of the setup (start of cycle)
i=1;                     %Value for plot
T_0 = 293;              %Starting temperature in [K], assumed for all fluids, gasses and materials
T_water = 293;              %Starting temperature water in [K]
T_al = 293;                 %Starting temperature aluminium in [K]
T_cu = 293;                 %Starting temperature copper in [K]
T_pvc = 293;                %Starting temperature PVC in [K]
T_pur = 293;                %Starting temperature polyurethane in [K]
T_air = 293;                %Starting temperature air in [K]
t = 0;                  %Time at the start [s]
t_final=12;                 %Time at end cycle [s]
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
hsv_pvc_d = 0.0022;     %Thickness of the HSV PCV [m]
r_outer_pvc = 0.055;    %Outer radius PVC tube [m]
r_inner_pvc = r_outer_pvc - hsv_pvc_d ;   %Inner radius PVC tube [m]
r_outer_double_pvc = r_outer_pvc + hsv_pvc_d; %Outer Radius of double-layered PVC tube (near endcaps) [m]
length_pvc = 0.130;     %Length/Height PVC tube [m] (this was shown to be slightly longer during first building phase with no real means to cut it down to 0.125m)
length_cap_pvc = 0.034; %Length/Height of a PVC Endcap [m]
epsilon_pvc = 0.9;      %Emissivity PVC [-]
rho_pvc = 1330;         %Density PVC [kg/m^3]
k_pvc = 0.19;           %Thermal conductivity PVC [W/(m K)]

%Polyurethane tubing
r_inner_pur = 0.004;    %Inner radius polyurethane tube [m]
r_outer_pur = 0.006;    %Outer radius polyurathane tube [m]
d_pur = 0.002;          %Wall Thickness polyurthane tube [m]
length_pur = 6;         %Length polyurethane tube [m]
epsilon_pur = 0.9;      %Emissivity polyurethane [-]
k_pur = 0.13;           %Thermal conducitivity polyurethane [W/(m K)]

%Aluminium plate (tape)
A_al = 1.085;           %Area aluminium plate [m^2]
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
k_air = 0.5;             %Placeholder

%General Setup Wood Properties
k_wood = 0.12;       %Thermal Conductivity of Pinewood [W/m k]                       
d_wood = 0.054;      %Wood border thcikness [m], simplification, in reality 2 mm air pocket within 
A_wood = 0.3107;     %Inner Surface Area wood layer [m], slight simplification

%Kingspan Insulation Properties
ks_k = 0.022;     % Thermal Conductivity Kingspan insulation layer[W/(m*k)]
ks_d = 0.04;      % thickness of Kingspan insulation layer [m]
ks_spacer_d = 0.01; % Thickness of a single Kingspan spacer layer [m]

%Other
sigma = 5.67*10^-8;     %Stefan Boltzmann constant [W/(m^2 K^4)]
flowrate = 3/60000;     %Flowrate pump [m^3/s]
E = 1000;               %Irradiance artificial sun [W/m^2]

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
% The following values, once again lifted from the earlier model, aim at extracting helpfull values for calculation of the resistive properties of air in between the HSV and its insulation.
% These values are simplified estimations.
    hsv_air_avgdist1 = ((0.010+0.00795)/2)*((length_pvc - (pvc_ec_inner_length * 2))/length_pvc) + ((0.0078+0.00575)/2)*((pvc_ec_inner_length * 2)/length_pvc); %Average distance of the heat reflector to the outer wall of the PVC body [m], USES PRE CALCULATED VALUES.  Accounts for the double layer PVC near the endcap, [m]
    hsv_air_avg_a1 = length_pvc * (hsv_al_width - (2*(sind(18)*0.005)));                                                                                                                %Surface area at the average point between pvc and reflector [m^2], (TODO implement avgdist 1)
    hsv_air_avgdist2 = 0.00905;                                                                                                                                                             %Average distance between reflector and kingspan therma plate [m]
    hsv_air_avg_a2 = length_pvc * (hsv_al_width + (2*(sind(18)*0.005)));                                                                                                                %Surface are at the average point between reflector and kingspan therma [m^2], (TODO implement avgdist 2)
    
%Aluminum Reflector Tape Properties
    d_al_t = 0.0001;                            % Thickness of an aluminum reflector segment [m]
    a_al_t = hsv_al_width * length_pvc;         % Area of a single heat reflector element. 10 are present in total [m^2]
    
% Convective heat transfer coefficients
h_air = 10;         %CVTH of still air [W/(m^2 K)], PLACEHOLDER
h_water = 0.5;      %CVTH of water [W/(m^2 K)], PLACEHOLDER
h_pvc = 0.5;        %CVTH PVC [W/(m^2 K)], PLACHEOLDER
h_al = 0.05;        %CVTH aluminium [W/(m^2 K)]
h_Kingspan = 0.5;   %CTVH of Kingspan-Therma insulation [W/(m^2 K)], PLACEHOLDER
h_pur = 0.5;        %CVTH of Polyurethane [W/(m^2 K)], PLACEHOLDER            

%Radiative heat transfer coefficients
h_r_pvc = 0.5;       %RDTH of pvc [W/(m^2 K)], PLACEHOLDER
h_r_al = 0.08;       %RDTH of aluminium-foil [W/(m^2 K)], PLACEHOLDER (but value is uncertain within a certain range) 
h_r_Kingspan = 0.5;  %RDTH of Kingspan-Therma insulation [W/(m^2 K)], PLACEHOLDER
h_r_pur = 0.5;       %RDTH of Polyurethane [W/(m^2 K)], PLACEHOLDER

%Overall Heat Transfer Coefficients, (may change due to change in CVTH)
U_cu = 1700;         %OHTC water in copper [W/(m^2 K)]
U_pur = 120.76;      %OHTC water in polyurethane [W/(m^2 K)]
U_pvc = 4.21;       %OHTC water in PVC [W/(m^2 K)]

%Other
V_system = V_cu+V_pvc+V_pur;                  %Volume of system  



%Geometry-based Therma Resistance
%check for reciprocal correctness
    %Radial
    R_a_cd = (log(r_outer_pvc / r_inner_pvc))/(2*pi * length_pvc * k_pvc); %Conductive transfer through pvc
    R_a_cv = 1/(h_pvc* (hsv_pvc_a_single+hsv_pvc_a_double)); %Convection into air from pvc
    R_b_cd = hsv_air_avgdist1 / (k_air * (10*hsv_air_avg_a1));  %Conduction through first still air pocket. simplified, assuming no convection.
    R_c_r = (10*a_al_t) * (h_r_al * h_al); %Radiation transfer into aluminum reflector CHECK THIS
    
    R_d_cd = d_al_t / (k_al * (10 * a_al_t) * length_pvc); %placeholder. Conduction through aluminum reflector. remember to account for the 10 panels
    R_d_cv = 1/((10*a_al_t) * h_al); %Convection from Aluminum plate into second air pocket
    R_d_r = R_c_r; %Radiation into second air pocket
    
    R_e_cd = hsv_air_avgdist2 / (k_air * (10*hsv_air_avg_a2)); %Conduction through second still air pocket
    R_e_cv = (10*hsv_therma_a_in) * (h_Kingspan + h_r_Kingspan);  %Convection into kingspan
    
    R_f_cd = ks_d * (ks_k * hsv_therma_a_avg * 10);  %Conduction through kingspan
    R_f_cv = (10*hsv_therma_a_out) * (h_r_Kingspan + h_Kingspan);  %Convection to outside
    
    
    %Endcap
    R_hsv_endcaps_1 = 0; %Conduction through endcap pvc material
    R_hsv_endcaps_2 = 0; %Convection into air pocket
    R_hsv_endcaps_3 = 0; %Conduction through air
    R_hsv_endcaps_4 = 0; %Convection into Kingspan therma
    R_hsv_endcaps_5 = 1 / (2 / (ks_d / (ks_k * (pvc_ec_a)))); %Conduction through kingspan therma
    R_hsv_endcaps_6 = 0; %Convection to outside
    
    %Total
        %Unsure as to whether or not to approach air pocket resistance as
        %series or parallel
    R_hsv_radial = R_a_cd + R_a_cv + R_b_cd + R_c_r + R_d_cd + R_d_cv + R_d_r + R_e_cd + R_e_cv + R_f_cd + R_f_cv;
    R_hsv_endcaps =  R_hsv_endcaps_1 + R_hsv_endcaps_2 + R_hsv_endcaps_3 + R_hsv_endcaps_4 + R_hsv_endcaps_5 + R_hsv_endcaps_6;
    %R_hsv = 1/(1/R_hsv_radial + 1/R_hsv_endcaps);
    R_hsv = 33.03; %substitute thermal resistance from 24-05 version.

%% Plotting info

while t<t_final
    %Water
    rho_water = (999.83953 + 16.945176 * (1.00024*T_water) - 7.9870401*10^-3 * (1.00024*T_water)^3 - 46.17046*10^-6* (1.00024*T_water)^3 +105.56302*10^-9 * (1.00024*T_water)^4 - 280.54253*10^-12 * (1.00024*T_water)^5)/(1+16.897850*10^-3 * (1.00024*T_water));   %Density water
    M_water = V_system*rho_water_20C;       %Volume of water inside system
    
    %Solar collector
    Q_rad_cu = E*length_cu*(r_outer_cu*2)*epsilon_paint;            %Heat addition radiation on copper tube [W]
    Q_rad_al = E*A_al*epsilon_paint;                                %Heat addition radiation on aluminium plate [W]
   
    Q_loss_conv_al = h_air*A_al*(T_al-T_air);                      %Heat loss convection aluminium plate [W]
    Q_loss_cond_al_cu = k_cu* A_exposed_cu*(T_al-T_cu)/(r_outer_cu*2-r_inner_cu*2);  %dit klopt niet hlml maar weten A tussen plaat en buis niet %Heat loss conduction aluminium plate [W]
    
    R_al = d_al/(k_al * A_al);                    %Thermal resistance aluminium plate 
    Q_loss_cond_al = (T_al-T_air)/R_al ;                                                              %Heat loss convection aluminium plate (other direction)
                            
    R_sol_air = 0.1;                              %Conductive thermal resistance air, PLACEHOLDER 
    R_sol_wood = d_wood/(k_wood * A_wood) ;       %Conductive thermal resistance wood setup
    
    Q_loss_rad_cu = sigma * epsilon_paint * A_outer_cu* (T_cu^4-T_air^4);    %Heat loss radiation copper tube [W]
    Q_loss_rad_al = sigma * epsilon_paint * A_al * (T_al^4-T_air^4);          %Heat loss radiation aluminium plate [W]
    
    T_al=T_al+(Q_rad_al-Q_loss_conv_al-Q_loss_cond_al-Q_loss_rad_al)/(M_al*c_al);  %Temperature of the aluminium plate [K]
    T_cu=T_cu+(Q_rad_cu-Q_loss_rad_cu-Q_loss_cond_al_cu)/(M_cu*c_cu);              %Temperature of the copper tube
    
    Q_conv_cu = U_cu * A_outer_cu * (T_cu - T_water);        %Convection copper-water
    Q_cu(i) = Q_conv_cu;
  
    
    %General Setup, (Adaptation of Heat Loss Stages in HSV)
    R_a_gs = 1/(A_inner_pur*h_water);                                         %Heat Absorption of Polyurethane
    R_b_gs = 1/((k_pur*A_inner_pur)/d_pur);                                   %Conduction through polyurethane tube
    R_c_gs = (r_outer_pur * 2 * pi * length_pur) * (h_r_pur * h_pur);         %Convection and Radiation into Outside air
    R_gs_total = R_a_gs + R_b_gs + R_c_gs;                                    %Total thermal resistance connection
   
    %Heat flux polyurethane tubes
    Q_gs = (T_water-T_0)/R_gs_total;
    
    %Heat flux HSV
        %if Q from other parts is known, t1 and t2 could be determined? (?  )
    
        
        
        
    %Temperature water
    T_water = T_0+((Q_conv_cu - Q_gs)/M_water*c_water);      %Final temperature water [K]
    y(i) = T_water;
    
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




   
