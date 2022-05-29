clear all, close all, clc

%% Notes for Version: 24 May, Revision 1
%This system simplifies the entie Solar-Thermal system to be a massive storage vessel with a certain heat input, and a certain heat loss.
%Time increments are 1.0 seconds and the simulation runs for the default 20 mins (1200 sec).
%Overall thermal model is COMPLETE but uses certain placeholders that can be easily calculated.
%Made for R2021a

%TODO
%Certain heat transfer coefficients still need implementation
%Heat Input currently uses a placeholder based on the exposed surface area of the copper tube (~0.1282 m^2)
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
flow_rate = 3 / 60;  %Volumetric flow of liquid pump, range 0.1/60 - 3.0/60, [L/s]
flow_rate_m3 = flow_rate / 1000; %Volumetric flow rate of pump, [m^3/s]
t_step = 1; % dt, time step interval [s]
sys_n_water_divisions = 200; %Self explanatory. Deprecated and no longer used in 13_05 onwards
sys_volume = 1.4; % Volume of entire system, [L]
sys_volume_m3 = sys_volume / 1000; %Volume of entire system, [m^3]
sys_mass = (sys_volume * 997)/1000; %Liquid mass based on volume, [kg]
sys_mass_sec = (flow_rate * 997)/1000; %Liquid mass based on volume, [kg]

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
c_cop = 385; %Specific Heat of copper, [J/Kg*k]

    %Thermal Transfer Coefficients 
    hsv_pvc_k = 0.19; %Thermal conductivity coefficient, for pvc should be in the 0.12 to 0.25 range, [W/(m*k)] 
    hsv_pvc_h = 5;  %Convection coefficient pvc into still air, [W/(m^2*k)]
    hsv_pvc_h_r = 0.5;                                          %placeholder
    hsv_air_k = 0.5;                                            %placeholder
    hsv_air_h = 0.5;                                            %placeholder
    hsv_air_h_r = 0.5;                                          %placeholder
    hsv_al_k = 237; %[W/(m*k)]
    hsv_al_h = 0.5;                                             %placeholder
    hsv_al_h_r = 0.08;                                          %placeholder but seems to be that this is another thing where it's a range of potential values. Anywhere from 0.04 to 0.1 should be reasonable in terms of the emissivity coefficient.
    hsv_therma_k = 0.022; %[W/(m*k)]
    hsv_therma_h = 0.5;                                         %placeholder
    hsv_therma_h_r = 0.5;                                       %placeholder
    hsv_water_k = 0.598; %[W/(m*k)]
    hsv_water_h = 0.5;                                          %placeholder
    col_cop_k = 400; %Thermal conductivity coefficient, [W/(m*k)]
    col_cop_h = 13.14; %Convection coefficient, [W/(m^2*k)]     


%------------------------------------------------
% Copper Solar Collector
    col_cop_reflectance = 0.05; %Reflectance Value on the surface of the Copper Tubing, dimensionless. Adjust between 0 and 1
    col_cop_ro = 0.006; %Outer radius, [m]
    col_cop_ri = 0.004; %Inner radius, [m]
    col_cop_e_roughness = 0.0015; %Pipe roughness e, [mm]
    col_cop_length = 6.85; %Length of exposed copper tube segments, [m]
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
    hsv_pvc_endcap_a = (hsv_pvc_ro)^2 * pi; %Surface area of a single endcap, [m^2]

    hsv_pvc_inner_endcap_length = hsv_pvc_endcap_length - hsv_pvc_thickness; %Length of the inner circular wall on the endcap. Used for thermal resistance calculation.
    hsv_pvc_d_single = hsv_pvc_thickness; %Used for R_b
    hsv_pvc_d_double = hsv_pvc_thickness * 2; %Used for R_b
    hsv_pvc_a_single = (hsv_pvc_length - (hsv_pvc_inner_endcap_length * 2)) * ((hsv_pvc_ri+hsv_pvc_ro)/2) * 2 * pi;
    hsv_pvc_a_double = (hsv_pvc_inner_endcap_length * 2) * ((hsv_pvc_ri+hsv_pvc_ro2)/2) * 2 * pi;
    hsv_al_width = 0.04017; % Width of Aluminum Tape reflector section, [m]
    hsv_al_d = 0.0001; %Thickness of the aluminum reflector segment, [m]
    hsv_al_a = hsv_al_width * hsv_pvc_length; %Area of a single heat reflector element. 10 are present in total, [m^2]
    hsv_therma_d = 0.04; %thickness of the kingspan therma insulation layer, [m]
%% The group of Values below are still somewhat questionable
    hsv_air_avgdist1 = ((0.010+0.00795)/2)*((hsv_pvc_length - (hsv_pvc_inner_endcap_length * 2))/hsv_pvc_length) + ((0.0078+0.00575)/2)*((hsv_pvc_inner_endcap_length * 2)/hsv_pvc_length); %USES PRE CALCULATED VALUES. Average distance of the heat reflector to the outer wall of the PVC body. Accounts for the double layer PVC near the endcap, [m]
    hsv_air_avg_a1 = hsv_pvc_length * (hsv_al_width - (2*(sind(18)*0.005))); %Surface area at the average point between pvc and reflector, ['m^2'] TODO implement avgdist 1
    hsv_air_avgdist2 = 0.00905; %Average distance between reflector and kingspan therma plate, [m]
    hsv_air_avg_a2 = hsv_pvc_length * (hsv_al_width + (2*(sind(18)*0.005)));  %Surface are at the average point between reflector and kingspan therma, [m^2] TODO implement avgdist 2
%%
    hsv_therma_width_in = (hsv_al_width + (2*(sind(18)*0.010))); %'width' of kingspan therma plate on the inside, [m]
    hsv_therma_width_out = (hsv_al_width + (2*(sind(18)*0.010))) + 2*(tand(18)*hsv_therma_d); %'width' of kingspan therma plate on the outside, [m]
    hsv_therma_a_in = hsv_pvc_length * hsv_therma_width_in; %Area of kingspan therma plate on inside, [m^2]
    hsv_therma_a_out = hsv_pvc_length * (hsv_therma_width_out); %Area of kingspan therma plate on outside, [m^2]
    hsv_therma_a_avg = (hsv_therma_a_out + hsv_therma_a_in)/2; %Used to calculate R_j, surface area in middle of therma plate, [m^2]




    
    %See Diagram in SSA 4 (Martijn) for schematic layout of resistances.
    %occurances of '10' are to account for decagon shape
    R_a = 1 / ((hsv_pvc_a_single + hsv_pvc_a_double) * hsv_water_h); %Absorption of heat energy from liq to pvd through convection
    R_b = 1 / (((hsv_pvc_k*hsv_pvc_a_single) / hsv_pvc_d_single) + ((hsv_pvc_k*hsv_pvc_a_double) / hsv_pvc_d_double)); %Complicated conduction process through PVC
    R_c = ((hsv_pvc_ro * 2 * pi) * hsv_pvc_length) * (hsv_pvc_h_r * hsv_pvc_h); %Emission through convection and radiation from pvc into still air
    R_d = hsv_air_avgdist1 / (hsv_air_k * (hsv_air_avg_a1)); %Conduction through still air
    R_e = hsv_al_a * (hsv_al_h_r * hsv_al_h); %Absorption through convection and Radiation into Al-foil
    R_f = hsv_al_d / (hsv_al_k * (hsv_pvc_ro * 2 * pi) * hsv_pvc_length); %Conduction through thin Al-foil
    R_g = R_e; %Same as R_g, but emission of heat energy from warm Al-foil through convection and radiation. Value unchanged.
    R_h = hsv_air_avgdist2 / (hsv_air_k * hsv_air_avg_a2 * 10); %Conduction through still air
    R_i = 10*(hsv_therma_a_in) * (hsv_therma_h + hsv_therma_h_r); %Convective and Radiation absorbtion by Kingspan therma panel
    R_j = hsv_therma_d * (hsv_therma_k * hsv_therma_a_avg *10 ); %Conduction through kingspan therma
    R_k = 10*(hsv_therma_a_out) * (hsv_therma_h_r + hsv_therma_h); %Convection and Radiation to Outside
    R_hsv_endcaps = 1 / (2 / (hsv_therma_d / (hsv_therma_k * (hsv_pvc_endcap_a)))); %Already accounts for the parallel resistance of both endcaps. TODO Might need to add the 2.2mm of PVC to this. Thermal Resistance HSV   
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


%------------------------------------------------
    %Mass Flow and Related
    col_cop_volume_interior = (col_cop_ri*col_cop_ri) * pi * col_cop_length; %Interior volume of the solar collector section. unused at the moment
%% Time-Dependent function 
for i = 1:t_final
    dQdt_in = col_cop_a_effective * P;
    Q_in(i) = dQdt_in*i; %For every entry. Delta-U since t=0 
    DBG_Q_in(i) = (sys_mass_sec/sys_mass)*(Q_in(i) / (sys_mass_sec*c_water)); %debug, just checking to see delta-T per timestep
    T_sys(i) = T_sys(i) + ((sys_mass_sec / sys_mass) * (Q_in(i) / (sys_mass_sec*c_water))); %System Temperature after heat energy increase
    
    dQdt_out(i) = ((T_sys(i)-T_amb)/R_hsv); %Essentially Heat flow. Add any extra terms for e.g. Feed and Return pipes here. Currently only HSV Present. This value should increase.
    Q_out_stor = Q_out_stor + dQdt_out(i); %Helper value, help with adding of total heat lost up until point i.
    Q_out(i) = Q_out_stor + dQdt_out(i); %For every entry. Delta-U since t=0 
    DBG_Q_out(i) = Q_out(i)/(sys_mass * c_water); %debug, just checking to see delta-T per timestep
    T_sys(i) = T_sys(i) - (Q_out(i) / (sys_mass*c_water)); %System Temperature after heat energy decrease


    %WIP section made by Stijn for heat loss in feed/return pipe sections
    Q_out_feed(i) = 2 * pi * feed_len * (T_sys(i) - T_amb) / (log(feedreturn_ro / feedreturn_ri) / feedreturn_k); %conductive heat loss through the pipe
    T_sys(i) = T_sys(i) - (Q_out_feed(i) / (sys_mass*c_water));
    
    Q_out_return(i) = 2 * pi * return_len * (T_sys(i) - T_amb) / (log(feedreturn_ro / feedreturn_ri) / feedreturn_k); %conductive heat loss through the pipe
    T_sys(i) = T_sys(i) - (Q_out_return(i) / (sys_mass*c_water));
end


%% Post Processing

%TODO implement multi-graph view of temperature and heat flow development over time.
plot(T_sys);
xlabel("Time Since Start");
ylabel("Liquid Temperature (k)");
