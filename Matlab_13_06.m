clear all, close all, clc

%% Notes for Version: 18 May, Revision 1
%This system simplifies the entie Solar-Thermal system to be a massive storage vessel with a certain heat input, and a certain heat loss.
%Time increments are 1.0 seconds and the simulation runs for the default 20 mins (1200 sec).
%Overall thermal model is COMPLETE but uses certain placeholders that can be easily calculated.
%Made for R2021a

%TODO
%Certain heat transfer coefficients still need implementation
%Heat Input currently uses a placeholder based on the exposed surface area of the copper tube (~0.1182 m^2)
%Feed and Return tubes need thermal resistance implementation
%Post Processing operations to create proper visualization.
%% Setup
% given starting variables
t_final = 1200; %Total amount of time for simulation, [s]
T_amb = 293; %Ambient Temperature, [k]
T_0 = 293; %Starting Temperature, equal to T_amb at start, [k]
T_sun = 333; %Assumes that the artificial sun will immediately reach 60 degrees celsius, [k]
P = 1000; %Power output of artificial sun, [W/m^2]

% System Variables
flow_speed = 3 / 60;  %Volumetric flow of liquid pump, range 0.1/60 - 3.0/60, [L/s]
t_step = 1; % dt, time step interval [s]
sys_n_water_divisions = 200; %Self explanatory. Deprecated and no longer used in 13_05 onwards
sys_volume = 1.4; % Volume of entire system, [L]
sys_mass = (sys_volume * 1000)/1000; %Liquid mass based on volume
sys_volume_water_divisions = sys_volume/sys_n_water_divisions; %Volume of each division. Deprecated and no longer used in 13_5 onwards.

% Temperature points. You can reference these later on for performance measuring.
dQdt_in = repmat(0, 1,t_final); %Create one-dimensional matrix to store heat flow (gain heat energy) at given i in.
Q_in = repmat(0, 1,t_final); %Create one-dimensional matrix to store total gained energy in.
T_sys = repmat(T_0, 1,t_final); %Create one-dimensional matrix to store liquid temperature in.
dQdt_out = repmat(0, 1,t_final); %Create one-dimensional matrix to store heat flow (lose heat energy) at given i in.
Q_out = repmat(0, 1,t_final); %Create one-dimensional matrix to store total lost energy in.
Q_out_stor = 0; %Declaration of helper value used later. See for-loop for implementation.
%% Material Properties
c_water = 4184; %Specific heat of water [J/kg*k]

%------------------------------------------------
% Copper Solar Collector
    col_cop_reflectance = 0.05; %Reflectance Value on the surface of the Copper Tubing, dimensionless. Adjust between 0 and 1.
    col_cop_k = 400; %Thermal conductivity coefficient, [W/(m*k)]
    col_cop_h = 13.14; %Convection coefficient, [W/(m^2*k)]
    col_cop_ro = 0.006; %Outer radius, [m]
    col_cop_ri = 0.004; %Inner radius, [m]
    col_cop_e_roughness = 0.0015; %Pipe roughness e, [mm]
    col_cop_length = 6.6; %Length of exposed copper tube segments, [m]

    c_cop = 385; %Specific Heat of copper, [J/Kg*k]
    col_cop_a = col_cop_ro*2*pi*col_cop_length; %Calculates the true surface area
    col_cop_a_exposed = col_cop_ro*pi* col_cop_length; %Calculates the surface area of the top half of the copper tubing, [m]
    col_cop_a_effective = col_cop_a_exposed * (1-col_cop_reflectance); %Applies reflectance parameter, [m]
    col_cop_mass = (col_cop_ro*2*pi- col_cop_ri*2*pi)*col_cop_length; %Mass of the copper tubing, [kg]
    
%Thermal Resistance Copper Solar Collector
    

%------------------------------------------------    
% Heat Storage Vessel
    hsv_pvc_thickness = 0.0022; %Wall Thickness of PVC, [m]
    hsv_pvc_length = 0.125; %Length of HSV, [m]
    hsv_pvc_ro = 0.055; %Outer Radius of HSV, [m]
    hsv_pvc_ri = 0.0528; %Inner Radius of HSV, [m]
    hsv_pvc_ro2 = 0.0572; %Outer Radius of double-layered parts near endcaps, [m]
    hsv_pvc_endcap_length = 0.034; %Length of endcap, [m]
    
    hsv_pvc_inner_endcap_length = hsv_pvc_endcap_length - hsv_pvc_thickness; %Length of the inner circular wall on the endcap. Used for thermal resistance calculation.
    hsv_pvc_d_single = hsv_pvc_thickness; %Used for R_b
    hsv_pvc_d_double = hsv_pvc_thickness * 2; %Used for R_b
    hsv_pvc_a_single = (hsv_pvc_length - (hsv_pvc_inner_endcap_length * 2)) * ((hsv_pvc_ri+hsv_pvc_ro)/2) * 2 * pi;
    hsv_pvc_a_double = (hsv_pvc_inner_endcap_length * 2) * ((hsv_pvc_ri+hsv_pvc_ro2)/2) * 2 * pi;
    
    hsv_al_a = 0.031 * 0.125; %Area of a single heat reflector element. 10 are present in total, [m^2]
    hsv_air_avgdist1 = ((0.010+0.00795)/2)*((hsv_pvc_length - (hsv_pvc_inner_endcap_length * 2))/hsv_pvc_length) + ((0.0078+0.00575)/2)*((hsv_pvc_inner_endcap_length * 2)/hsv_pvc_length); %Average distance of the heat reflector to the outer wall of the PVC body. Accounts for the double layer PVC near the endcap, [m]
    hsv_air_avgdist2 = 0;
    
    hsv_pvc_k = 0.19; %Thermal conductivity coefficient, for pvc should be in the 0.12 to 0.25 range, [W/(m*k)] 
    hsv_pvc_h = 5;  %Convection coefficient pvc into still air, [W/(m^2*k)]
    hsv_pvc_h_r = 0.5; %placeholder
    hsv_air_k = 0.5; %placeholder
    hsv_air_h = 0.5; %placeholder
    hsv_air_h_r = 0.5; %placeholder
    hsv_al_k = 237; 
    hsv_al_h = 0.5; %placeholder
    hsv_al_h_r = 0.5; %placeholder
    hsv_therma_k = 0.022; 
    hsv_therma_h = 0.5; %placeholder
    hsv_therma_h_r = 0.5; %placeholder
    hsv_water_k = 0.598; 
    hsv_water_h = 0.5; %placeholder
    hsv_pvc_surf_dub = 0.021099; %surface area of radial part of HSV where the thickness of PVC is 4.4mm, [m^2]
    hsv_pvc_surf_sin = 0.020369; %surface area of radial part of HSV where the thickness of PVC is 2.2mm, [m^2]
    
    %See Diagram in SSA 4 (Martijn) for schematic layout of resistances
    R_a = 1 / ((hsv_pvc_a_single + hsv_pvc_a_double) * hsv_water_h); %Absorption of heat energy from liq to pvd through convection
    R_b = 1 / (((hsv_pvc_k*hsv_pvc_a_single) / hsv_pvc_d_single) + ((hsv_pvc_k*hsv_pvc_a_double) / hsv_pvc_d_double)); %Complicated conduction process through PVC
    R_c = ((hsv_pvc_ro * 2 * pi) * hsv_pvc_length) * (hsv_pvc_h_r * hsv_pvc_h); %Emission through convection and radiation from pvc into still air
    R_d = hsv_air_avgdist1 / (hsv_air_k * 0.048449); %Conduction through still air
    R_e = hsv_al_a * (hsv_al_h_r * hsv_al_h); %Absorption through convection and Radiation into Al-foil
    R_f = 0.0001 / (hsv_al_k * 0.03875); %Conduction through thin Al-foil
    R_g = R_e; %Same as R_g, but emission of heat energy from warm Al-foil through convection and radiation. Value unchanged.
    R_h = 0.00905 / (hsv_air_k * 0.0555); %Conduction through still air
    R_i = 0.05906* (hsv_therma_h + hsv_therma_h_r); %Convective and Radiation absorbtion by Kingspan therma panel
    R_j = 0.04 * (hsv_therma_k * 0.0747699); %Conduction through kingspan therma
    R_k = 0.090477 * (hsv_therma_h_r + hsv_therma_h); %Convection and Radiation to Outside
    R_hsv_endcaps = 1 / (2 / (0.04 / (hsv_therma_k * 0.0095033))); %Already accounts for the parallel resistance of both endcaps. TODO Might need to add the 2.2mm of PVC to this.
%Thermal Resistance HSV   
R_hsv_rad = R_a + R_b + R_c + R_d + R_e + R_f + R_g + R_h + R_i + R_j + R_k; %summate all elements of the radial section together.
R_hsv = 1/(1/R_hsv_rad + 1/R_hsv_endcaps); %Add radial and endcap resistances as parallel resistance set.


%------------------------------------------------
    %Geometric Properties of Feed and Return pipe
    feedreturn_ri = 0.004; %inner radius of feed/return tube
    feedreturn_ro = 0.006; %outer radius of feed/return tube
    feed_len = 3; %length of feed tube
    return_len = 3; %length of return tube
    feedreturn_k = 0.13; %Thermal Conductivity coeff of the feed and return tubes, [W/(m*k)]
%Thermal Resistance Feed Pipe


%Thermal Resistance Return Pipe

Q_test = 10;

%% Time-Dependent function 
%for i = 1:t_final
    %dQdt_in = 100; %Fixed placeholder input power heat flow in watt. With real values this should decrease slightly over time as heat flux over copper tube decreases with increasing liq temp. TODO improve
   %Q_in(i) = dQdt_in*i; %For every entry. Delta-U since t=0 
   % DBG_Q_in(i) = Q_in(i) / (sys_mass*c_water); %debug, just checking to see delta-T per timestep
    %T_sys(i) = T_sys(i) + (Q_in(i) / (sys_mass*c_water)); %System Temperature after heat energy increase
    
   % dQdt_out(i) = ((T_sys(i)-T_amb)/R_hsv); %Essentially Heat flow. Add any extra terms for e.g. Feed and Return pipes here. Currently only HSV Present. This value should increase.
   % Q_out_stor = Q_out_stor + dQdt_out(i); %Helper value, help with adding of total heat lost up until point i.
   % Q_out(i) = Q_out_stor + dQdt_out(i); %For every entry. Delta-U since t=0 
   % DBG_Q_out(i) = Q_out(i)/(sys_mass * c_water); %debug, just checking to see delta-T per timestep
   % T_sys(i) = T_sys(i) - (Q_out(i) / (sys_mass*c_water)); %System Temperature after heat energy decrease
    
%end

   
    
 

t = linspace(0,10,10);
y0 = 0.0;
u = 1.0;
Temperature = ode23(@(t,y)derrivative(t,y,u),t,y0);
Full_temperature = Temperature.x;
y = Temperature.y;
plot(Full_temperature,y)


    function dydt = derrivative(t,y,u)
 radiation = 5;
 Q = 2.0;
 
 dydt= (-y + Q*u/radiation);
 
 end 
   
    
%% Post Processing
%TODO implement multi-graph view of temperature and heat flow development
%over time.

