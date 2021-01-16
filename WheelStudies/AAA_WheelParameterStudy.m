% Wheel Parameter Study
% by Austin Lines
% Last updated: 2021-01-15
% Referenced in Chapter 3 of PhD Thesis

% Shear and normal stress of a rigid tire
% Equations from Wong

% Iterative determination of theta1 given slip ratio
% Computation of resistance and total torque

%%
clear
close all
clear global

global g Lf Lr hcg Rw mg m muf mur ms LfpLr Iwf Iwr e track_f tr ho Ixxs Ixzs Ixz Izz mhur mhuf ... 
    hfLrp hrLfm mshcg mshcgg mshchpho tfo2 tro2 Fzf_1 Fzr_1 beta gam ho6 ho8 ho2 ho3 Bw;

%% SETTINGS

% FORCE-SLIP CURVE
% 1 - ON
% 0 - OFF
slipCalc = 0;

% PRESSURE-SINKAGE MODEL
% 1 - Bekker model
% 2 - Preston-Thomas model
% 3 - Snow
% 4 - Reece
model = 2;

% SHEAR STRESS-SHEAR DISPLACEMENT MODEL
% 1 - Poorly Bonded (Wong - Theory of Ground Vehicles - pg. 144, A.)
% 2 - Well Bonded   (Wong - Theory of Ground Vehicles - pg. 149, C.)
shear_model = 1;

% PATCH
% 0 - max pressure at theta = 0, no constant pressure patch
% 1 - constant pressure contact patch
% 2 - Wong-Reece - maximum pressure somewhere in middle of wheel
patch = 2;
z_const = .056;         % meters

% SLIP-SINKAGE
% 0 - OFF
% 1 - Lyasko's model of (1 + i)/(1 - 0.5*i)
% 2 - Lyasko's model of z = z0 + i*S*z0
% 3 - Ding's Model
slip_sink = 2;
S = 0.6/0.33;

if slip_sink == 2
    % FLATTEN AT HIGH SINKAGE
    % 0 - OFF
    % 1 - Flatten
    flatten = 1;
end
    
% BULLDOZING RESISTANCE
% 0 - OFF
% 1 - ON
bulldozing = 1;

% ROBOT
% 1 - Default
% 2 - Cool Robot
% 3 - FrostyBoy
robot = 3;

% TERRAIN DAMPER
% 0 - OFF
% 1 - ON
damp_terrain = 1;
b_terrain = 1/200;
sink_max = .001;  %.001 standard

% THETA2 CONTROL
% 0 - OFF
% 1 - ON
% 2 - Ishigami equation
slip_theta2 = 0;
lambda = .11;        % Lambda for Ishigami equation

ss_nom = 0.1;   % nominal stepsize for integration of shear and normal stress

%% List of common terrain parameters
% To use a set of terrain parameters below, copy and paste the terrain
% parameters into the next section

run = false;
if run == true
    
    % Clayey soil
    phi = 13*pi/180;    % angle of internal shearing resistance
    c = 4.14;           % cohesion of soil
    Ks =0.01;           % shear deformation modulus
    kc = 13.19;
    kphi = 692.15;      % pressure sinkage parameters
    n = 0.5;
    c1 = 0.18;
    c2 = 0.32;          % these are probably not needed for deformable tire

    % Dry sand from Iagnemma et al. 2004
    phi = 30*pi/180;    % angle of internal shearing resistance
    c = 1.0;            % cohesion of soil  kPa 
    Ks =0.025;          % shear deformation modulus (m)
    kc = 0.9;           % kPa/m^n-1 
    kphi = 1523.4;      % pressure sinkage parameters  kPa/m*n
    n = 1.1;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Sandy loam  from Wong p. 136
    phi = 33.7*pi/180;  % angle of internal shearing resistance
    c = 3.3;            % cohesion of soil  kPa 
    Ks = 0.025;          % shear deformation modulus (m)
    kc = 74.6;          % kPa/m^n-1 
    kphi = 2080;        % pressure sinkage parameters  kPa/m*n
    n = 1.1;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Clayey soil
    phi = 13*pi/180;    % angle of internal shearing resistance
    c = 4.14;           % cohesion of soil
    Ks = 0.01;           % shear deformation modulus
    kc = 13.19;
    kphi = 692.15;      % pressure sinkage parameters
    n = 0.5;
    c1 = 0.18;
    c2 = 0.32;          % these are probably not needed for deformable tire

    % Snow from Wong p. 136 (Sweden)
    phi = 20.7*pi/180;  % angle of internal shearing resistance
    c = 6;              % cohesion of soil  kPa 
    Ks =0.02;           % shear deformation modulus (m)
    kc = 10.55;         % kPa/m^n-1 
    kphi = 66.08;       % pressure sinkage parameters  kPa/m*n
    n = 1.44;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Snow from Wong p.136 (Harrison)
    phi = 23.2*pi/180;  % angle of internal shearing resistance
    c = .62;            % cohesion of soil  kPa 
    Ks =0.02;           % shear deformation modulus (m)
    kc = 2.49;          % kPa/m^n-1 
    kphi = 245.9;       % pressure sinkage parameters  kPa/m*n
    n = 1.6;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Snow from Wong p.136 (US)
    phi = 19.7*pi/180;  % angle of internal shearing resistance
    c = 1.03;           % cohesion of soil  kPa 
    Ks = 0.02;          % shear deformation modulus (m)
    kc = 4.37;          % kPa/m^n-1 
    kphi = 196.72;      % pressure sinkage parameters  kPa/m*n
    n = 1.6;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Snow from Lever
    phi = 20.7*pi/180;  % angle of internal shearing resistance
    c = .5;             % cohesion of soil  kPa 
    Ks = 0.02;          % shear deformation modulus (m)
    kc = 90;            % kPa/m^n-1 
    kphi = 90;          % pressure sinkage parameters  kPa/m*n
    n = 1.6;
    c1 = 0.3;
    c2 = 0.4;           % these are probably not needed for deformable tire

    % Lean clay
    phi = 20*pi/180;    % angle of internal shearing resistance
    c = 68.95;          % cohesion of soil
    K = 0.0001;         % shear deformation modulus
    Ks = K;
    kc = 16.42;
    kphi = 1724.69;     % pressure sinkage parameters
    n = 0.2;
    c1 = 0.43;
    c2 = 0.32;

end
%% Define terrain parameters to be used in study

% BEKKER
phi_d = 9;
phi = phi_d*pi/180;     % angle of internal shearing resistance
c = 1.2;                % cohesion of snow  kPa
Ks = .03;               % shear deformation modulus (m)
                        % from Wong - "for fresh snow, values range from 2.5 to 5 cm" 
kc = 17.16;             % kPa/m^n-1 
kphi = 8.038;           % pressure sinkage parameters  kPa/m*n
n = 3.596;

% WONG, PRESTON-THOMAS
zm = .54;           % meters
km = 40;            % kPa/m
c1 = 0.18;          % Wong and Reece 1967, loose sand
c2 = 0.32;          % Wong and Reece 1967, loose sand

% SNOW
kp1 = 16.3;         % kPa
kp2 = 0;            % kPa/m
kz1 = .033;         % meters
kz2 = 0;            % meter^2

% REECE
% kc = 50;
% n = 2;
% Ks = 0.061;

% WONG, shear stress-shear displacement, Model C
Kw = 0.022;
Kr = 0.66;

% Soil bearing capacity factors
[Kc,Ky] = soil_factors(phi);
gamma = 1961;        % N/m^3

%% Vehicle parameters

if robot == 1   % Default
    g = 9.81;               % Gravity constant, m/s^2
    mv = 9.6071655;         % Vehicle mass, kg
    m = mv;
    Rw = 5/100;             % Wheel radius in m
    mw = 0.28;              % Wheel mass, kg
    Iw = mw*Rw^2;           % Wheel moment of inertia, kg m^2
    Iwf = Iw;               % Front wheel moment of interia, kg m^2
    Iwr = Iw;               % Rear wheel moment of interia, kg m^2
    muf = mw*2;             % Total front wheel mass, kg
    mur = mw*2;             % Total rear wheel mass, kg
    ms = m - muf - mur;     % Vehicle body mass (w/o wheels), kg
    tr = 0.43722;           % Track width, m
    track_f = tr;           % Front track width, m, non-trivial to implement as not =
    Lr = 0.184361;          % Distance from CG to rear axle centerline, m
    Lf = 0.216237;          % Distance from CG to front axle centerline
    LfpLr = Lr + Lf;        % Total vehicle length between axle centers, m
    hcg = 0.05;             % Height from wheel axis to vehicle CG
    hf = 0.05;              % Height of front roll center w.r.t. ground, m
    hr = 0.05;              % Height of rear roll center w.r.t. ground, m
    huf = Rw;               % Front unsprung mass height (wheel center) w.r.t. ground, m
    hur = Rw;               % Rear unsprung mass height (wheel center) w.r.t. ground, m
    ho = (hf*Lr + hr*Lf)/LfpLr;
    e = 0;                  % Horizontal distance from body CG to total vehicle CG
    Izz = .23687393;        % Z-axis (yaw) moment of inertia, kg*m^2
    mhur = mur*hur;
    mhuf = muf*huf;
    hfLrp = hf*(Lr + e)/LfpLr;
    hrLfm = hr*(Lf - e)/LfpLr;
    mshcg = ms*hcg;
    mshcgpho = ms*(hcg + ho);
    mshcgg = mshcg*g;
    tfo2 = track_f/2;
    tro2 = tr/2;
    mg = m*g;               % Total vehicle normal force Newtons
    Fz = mv*g;              % Total vehicle normal force Newtons
    
    Fzf_1 = Lr*m*9.81/(Lf + Lr);
    Fzr_1 = Lf*m*9.81/(Lf + Lr);
    Fzfl = (Lr*m*9.81/(Lf + Lr))/2; %Front left normal force, N
    Fzfr = (Lr*m*9.81/(Lf + Lr))/2; %Front right normal force, N
    Fzrl = (Lf*m*9.81/(Lf + Lr))/2; %Rear left normal force, N
    Fzrr = (Lf*m*9.81/(Lf + Lr))/2; %Rear right normal force, N
    
    Bw = 0.05; % damping coefficient for tires

elseif robot == 2 % Cool Robot
    g = 9.81;           % Gravity constant, m/s^2
    mv = 60;            % Vehicle mass, kg
    m = mv;
    Rw = .4572/2;       % Wheel radius in m (Cool Robot)
    b = 0.127;          % Wheel width (Cool Robot)
    mw = 11/2.2;        % Wheel mass, kg
    Iw = mw*Rw^2;       % Wheel moment of inertia, kg m^2
    Iwf = Iw;           % Front wheel moment of interia, kg m^2
    Iwr = Iw;           % Rear wheel moment of interia, kg m^2
    muf = mw*2;         % Total front wheel mass, kg
    mur = mw*2;         % Total rear wheel mass, kg
    ms = m - muf - mur; % Vehicle body mass (w/o wheels), kg
    tr = .9745;         % Track width, m
    track_f = tr;       % Front track width, m, non-trivial to implement as not =
    Lr = .3155;         % Distance from CG to rear axle centerline, m
    Lf = .3155;         % Distance from CG to front axle centerline, m
    LfpLr = Lr + Lf;    % Total vehicle length between axle centers, m
    hcg = 0.05;         % Height from wheel axis to vehicle CG
    hf = 0.05;          % Height of front roll center w.r.t. ground, m
    hr = 0.05;          % Height of rear roll center w.r.t. ground, m
    huf = Rw;           % Front unsprung mass height (wheel center) w.r.t. ground, m
    hur = Rw;           % Rear unsprung mass height (wheel center) w.r.t. ground, m
    ho = (hf*Lr + hr*Lf)/LfpLr;
    e = 0;              % Horizontal distance from body CG to total vehicle CG
    Izz = 11.4586;      % Z-axis (yaw) moment of inertia, kg*m^2
    mhur = mur*hur;
    mhuf = muf*huf;
    hfLrp = hf*(Lr + e)/LfpLr;
    hrLfm = hr*(Lf - e)/LfpLr;
    mshcg = ms*hcg;
    mshcgpho = ms*(hcg + ho);
    mshcgg = mshcg*g;
    tfo2 = track_f/2;
    tro2 = tr/2;
    mg = m*g;           % Total vehicle weight, N
    Fz = mv*g;          % Total vehicle normal force Newtons
    
    % Normal forces
    Fzf_1 = Lr*m*9.81/(Lf + Lr);
    Fzr_1 = Lf*m*9.81/(Lf + Lr);
    Fzfl = (Lr*m*9.81/(Lf + Lr))/2; %Front left normal force, N
    Fzfr = (Lr*m*9.81/(Lf + Lr))/2; %Front right normal force, N
    Fzrl = (Lf*m*9.81/(Lf + Lr))/2; %Rear left normal force, N
    Fzrr = (Lf*m*9.81/(Lf + Lr))/2; %Rear right normal force, N
    
    Bw = 0.0005;        % Damping coefficient for tires
    
elseif robot == 3       % FrostyBoy Robot
    g = 9.81;           % Gravity constant, m/s^2
    mv = 90;            % Vehicle mass, kg
    m = mv;
    Rw = .265;          % Wheel radius in m
    b = .45;            % Wheel width
    mw = 13/2.2;        % Wheel mass, kg
    Iw = mw*Rw^2;       % Wheel moment of inertia, kg m^2
    Iwf = Iw;           % Front wheel moment of interia, kg m^2
    Iwr = Iw;           % Rear wheel moment of interia, kg m^2
    muf = mw*2;         % Total front wheel mass, kg
    mur = mw*2;         % Total rear wheel mass, kg
    ms = m - muf - mur; % Vehicle body mass (w/o wheels), kg
    tr = 1.14224;       % Track width, m
    track_f = tr;       % Front track width, m, non-trivial to implement as not =
    wheelbase = .579+.515;
    Lr = wheelbase/2;   % Distance from CG to rear axle centerline, m
    Lf = wheelbase/2;   % Distance from CG to front axle centerline, m
    LfpLr = Lr + Lf;    % Total vehicle length between axle centers, m
    hcg = 0.140;        % Height from wheel axis to vehicle CG
    hf = 0.05;          % Height of front roll center w.r.t. ground, m
    hr = 0.05;          % Height of rear roll center w.r.t. ground, m
    huf = Rw;           % Front unsprung mass height (wheel center) w.r.t. ground, m
    hur = Rw;           % Rear unsprung mass height (wheel center) w.r.t. ground, m
    ho = (hf*Lr + hr*Lf)/LfpLr;
    e = 0;              % Horizontal distance from body CG to total vehicle CG
    Izz = 36.05;        % Z-axis (yaw) moment of inertia, kg*m^2
    mhur = mur*hur;
    mhuf = muf*huf;
    hfLrp = hf*(Lr + e)/LfpLr;
    hrLfm = hr*(Lf - e)/LfpLr;
    mshcg = ms*hcg;
    mshcgpho = ms*(hcg + ho);
    mshcgg = mshcg*g;
    tfo2 = track_f/2;
    tro2 = tr/2;
    mg = m*g;           % Total vehicle weight, N
    Fz = mv*g;          % Total vehicle normal force, N

    % Normal forces
    Fzf_1 = Lr*m*9.81/(Lf + Lr);
    Fzr_1 = Lf*m*9.81/(Lf + Lr);
    Fzfl = (Lr*m*9.81/(Lf + Lr))/2; %Front left normal force, N
    Fzfr = (Lr*m*9.81/(Lf + Lr))/2; %Front right normal force, N
    Fzrl = (Lf*m*9.81/(Lf + Lr))/2; %Rear left normal force, N
    Fzrr = (Lf*m*9.81/(Lf + Lr))/2; %Rear right normal force, N
    
    % Motor parameters
    Tw_rad = -10.0507;      % Nm/(rad/sec) at wheel (from 'RP34_MotorCalcs.m')
    V_ref = 48;             % reference voltage
    w_noload = 7000;        % RPM
    w_noload_rad = 18.3260; % no load speed at wheel [rad/sec]
    gearbox_efficiency = 0.8;
    T_max = 88.64;             % 88.64 = 2.77*40*gearbox_efficiency;     
    
    Bw = 0.3269;        % Damping coefficient for tires
    
end

%% Initialize parameters that are not varying

lt_nom = 0.05;
min_flag = 0;
theta1 = 0.18;
slip = 0.25;
W_est = 0;
W = Fz/4;

%% Pick two parameters to vary and uncomment them from the following list:

slip_range = linspace(0,1,3);                          % Slip range
% W_range = linspace(Fz/4,Fz/2,11);                       % Normal force
Rw_range = linspace(0.2,0.5,3);                         % Wheel radius
% b_range = linspace(0.1,.4,11);                          % Wheel width
% km_range = linspace(20,65,11);                          % Terrain stiffness (Preston-Thomas)
% zm_range = linspace(0.02,0.15,11);                      % Max sinkage (Preston-Thomas)
% kc_range = linspace(10,100,11);                         % Cohesive modulus (Bekker)
% kphi_range = linspace(10,150,11);                       % Friction modulus (Bekker)
% n_range = linspace(0.8,3,11);                           % Exponent (Bekker)
% phi_range = linspace(9,23,15);                          % Angle of internal friction
% Ks_range = linspace(0.025,0.05,11);                     % Shear deformation modulus (Poorly Bonded)
% theta2_range = linspace(deg2rad(-30),deg2rad(0),11);    % Wheel exit angle
% c_range = linspace(0.6,1.8,11);                         % Cohesion
% Kr_range = linspace(0.05,2,11);                         % Ratio of residual shear stress to maximum shear stress (Well Bonded)
% Kw_range = linspace(0.07,0.2,11);                       % Shear displacement resulting in max shear stress (Well Bonded)
% S_range = linspace(.5,.6/.33,11);                       % Slip-sinkage factor


%% Ensure that the two parameters that are not commented out above are
% replaced in all instances below (use "Find & Replace" function):

forces = zeros(length(slip_range),length(Rw_range),18);
for j = 1:length(slip_range)
    tic
%     W = W_range(j);
    slip = slip_range(j);
%     b = b_range(j);
%     disp(num2str(j));
%     Ks = Ks_range(j);
%     phi = deg2rad(phi_range(j));
for i = 1:length(Rw_range)
%     slip=slip_range(i);
%     W = W_range(i);
%     b = b_range(i);
    Rw = Rw_range(i);
%     km = km_range(i); %zm = -0.0054*km+0.6455;
%     zm = zm_range(i);
%     kc = kc_range(i);
%     kphi = kphi_range(i);
%     phi = deg2rad(phi_range(i));
%     Ks = Ks_range(i);
%     theta2 = theta2_range(i);
%     c = c_range(i);
%     Kr = Kr_range(i);
%     Kw = Kw_range(i);
%     n = n_range(i);
%     c1 = c1_range(i);
%     S = S_range(i);

    z_offset = 0;
    W_est = 0;
    stepsize = 0.01;
    [Kc,Ky] = soil_factors(phi);

while (W_est > (W + 0.5)/1000 )   || (W_est <  (W - 0.5)/1000)
if patch == 0
    thetam = 0;
    theta2 = 0;
elseif patch == 1
    thetam = real(acos((z_const/Rw)+cos(theta1)));
elseif patch == 2
    if slip>0
        thetam = (c1+c2*slip)*theta1;
    else
        thetam = c1*theta1;
    end
    if slip > .5 && (slip_theta2 == 1)
        theta2 = deg2rad(lambda*(slip-.5)-15);
    else
        theta2 = deg2rad(-15);
    end
end

% Find W by integrating normal stress and shear stress, break into two regions theta_1 to
% theta_m and theta_m to theta_2
W1 = integral(@(xt)normal_stress_lines(xt,theta1,thetam,kc,b,kphi,n,Rw,'normal',model,km,zm,kp1,kp2,kz1,kz2,z_offset),thetam,theta1);
W3 = integral(@(xt)shear_stress_lines(xt,theta1,thetam,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'normal',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),thetam,theta1);
if patch == 0
    W2 = 0;
    W4 = 0;
elseif patch == 1
    if model == 2
        W2 = km*zm*(-log(1-(z_const/zm)));
    else
        W2 = (kc + kphi*b)*((Rw/b).^n)*((1 - cos(theta1)).^n);
    end
elseif patch == 2
    W2 = integral(@(xt)normal_stress_lines_2(xt,theta1,thetam,theta2,kc,b,kphi,n,Rw,'normal',model,km,zm,kp1,kp2,kz1,kz2,z_offset),theta2,thetam);
    W4 = integral(@(xt)shear_stress_lines_2(xt,theta1,thetam,theta2,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'normal',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),theta2,thetam);
end

W_est = Rw*b*(W1+W2+W3+W4);
W_error = (W/1000 - W_est);
theta1 = theta1 + stepsize*W_error;
if theta1 >= (pi/2)
    theta1 = pi/2;
    z_offset = z_offset + stepsize*W_error/10;
end

end

% Compute sinkage
if slip_sink == 1
    Kss = (1+slip)/(1-(0.5*slip));
    z0 = (1 - cos(theta1))*Rw;
    sinkage_slip = Kss*z0;
    theta_slip = acos(1 - (sinkage_slip/Rw));
    theta1 = theta_slip;
elseif slip_sink == 2
    z0 = (1 - cos(theta1))*Rw;
    if slip>.3 && flatten ==1
        factor = ((S*.5)./(1+exp(-8*(slip-.25))));
    else
        factor = S*slip;
    end
    sinkage_slip = z0 + factor*z0;
    theta_slip = acos(1 - (sinkage_slip/Rw));
    theta1 = theta_slip;
elseif slip_sink == 3
    z0 = (1 - cos(theta1))*Rw;
    sinkage_slip = z0 - ding_slip_sink(thetam,theta2,phi,slip);
    theta_slip = acos(1 - (sinkage_slip/Rw));
    theta1 = theta_slip;
end

sinkage = (1 - cos(theta1))*Rw + z_offset;


% Compute resistance
R1 = integral(@(xt)normal_stress_lines(xt,theta1,thetam,kc,b,kphi,n,Rw,'resist',model,km,zm,kp1,kp2,kz1,kz2,z_offset),thetam,theta1);
R3 = integral(@(xt)shear_stress_lines(xt,theta1,thetam,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'resist',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),thetam,theta1);
if patch == 0 || patch == 1
    R2 = 0;
    R4 = 0;
elseif patch == 2
    R2 = integral(@(xt)normal_stress_lines_2(xt,theta1,thetam,theta2,kc,b,kphi,n,Rw,'resist',model,km,zm,kp1,kp2,kz1,kz2,z_offset),theta2,thetam);
    R4 = integral(@(xt)shear_stress_lines_2(xt,theta1,thetam,theta2,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'resist',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),theta2,thetam);
end

% Bulldozing Resistance
if bulldozing == 1
    z_b = sinkage - (1-cos(thetam))*Rw;
    Rb = b*(0.667*z_b*c*Kc + 0.5*(z_b^2)*gamma*Ky);   % Gee-Clough, pg 162
    DP = Rw*b*(R3 + R4 - R1 - R2)*1000 - (Rb);
    R_terrain = Rw*b*(R1+R2)*1000 + Rb;
else
    DP = Rw*b*(R3 + R4 - R1 - R2)*1000;
    R_terrain = Rw*b*(R1+R2)*1000;
    Rb = 0;
end

% Compute total torque
T1 = integral(@(xt)shear_stress_lines(xt,theta1,thetam,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'torque',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),thetam,theta1);
if patch == 0 || patch == 1
    T2 = 0;
else
    T2 = integral(@(xt)shear_stress_lines_2(xt,theta1,thetam,theta2,kc,b,kphi,n,Rw,phi,c,Ks,Kr,Kw,slip,'torque',model,km,zm,kp1,kp2,kz1,kz2,0,min_flag,shear_model,z_offset),theta2,thetam);
end

T = Rw*Rw*b*(T1 + T2)*1000;

% Tractive effort
Ft = Rw*b*(R3+R4)*1000;

%                   1   2  3  4  5  6 7  8  9  10    11     12 13 14 15     16  17 18    
forces(j,i,:) = [theta1 W1 W2 W3 W4 W R1 R2 R3 R4 R_terrain T1 T2 T  Ft sinkage DP Rb];
end
toc
end

%% Plot comparisons sets of data over the two independent variables' ranges
% Some example plots are shown below, but need to be modified according to
% which independent variables were chosen

% The matrix "forces" consists of the following datasets:
% (j,i,1) - front sinkage angle (theta_1)
% (j,i,2) - normal force supporting wheel in z (theta_1 to theta_m)
% (j,i,3) - normal force supporting wheel in z (theta_m to theta_2)
% (j,i,4) - shear force supporting wheel in z (theta_1 to theta_m)
% (j,i,5) - shear force supporting wheel in z (theta_m to theta_2)
% (j,i,6) - total force supporting wheel in z (theta_1 to theta_2)
% (j,i,7) - normal force resisting wheel in x (theta_1 to theta_m)
% (j,i,8) - normal force resisting wheel in x (theta_m to theta_2)
% (j,i,9) - shear force resisting wheel in x (theta_1 to theta_m)
% (j,i,10) - shear force resisting wheel in x (theta_m to theta_2)
% (j,i,11) - total force resisting wheel in x (theta_1 to theta_2 + bulldozing)
% (j,i,12) - resistance torque on wheel (theta_1 to theta_2)
% (j,i,13) - resistance torque on wheel (theta_m to theta_2)
% (j,i,14) - total resistance torque on wheel (theta_1 to theta_2)
% (j,i,15) - tractive effort
% (j,i,16) - sinkage
% (j,i,17) - drawbar pull
% (j,i,18) - bulldozing resistance


% Example plots:

% Terrain resistance - width vs radius
% figure()
% hold on
% % for o = 1:length(Rw_range)
% plot(b_range-b_range(1),forces(:,1,11))
% plot(Rw_range-Rw_range(1),forces(1,:,11))
% % end
% xlabel('Change in Wheel Dimension')
% ylabel('Terrain Resistance [N]')
% legend('\Delta b','\Delta R');
% set(gcf,'Position',[675 40 560 420])

% Resistive torque - radius
% figure()
% hold on
% % plot(forces_b(1,:,16),forces_b(1,:,14)*1000,'+');
% % plot(forces_b(5,:,16),forces_b(5,:,14)*1000,'x');
% % plot(forces_b(10,:,16),forces_b(10,:,14)*1000,'*');
% plot(forces(1,:,16),forces(1,:,14),'o');
% plot(forces(2,:,16),forces(2,:,14),'s');
% plot(forces(3,:,16),forces(3,:,14),'d');
% guide{1} = '\Delta b, W = 60 kg';
% guide{2} = '\Delta b, W = 90 kg';
% guide{3} = '\Delta b, W = 120 kg';
% legend(guide);
% xlabel('Sinkage [m]');
% ylabel('Resistance Torque [N-m]');
% title({'Torque vs Sinkage','Sensitivity to Wheel Width and Diameter'});
% set(gcf,'Position',[675 40 560 420])
% set(gca,'XDir','reverse')

% figure()
% hold on
% plot(forces(1,:,16),forces(1,:,17)./forces(1,:,14),'o');
% plot(forces(2,:,16),forces(2,:,17)./forces(2,:,14),'s');
% plot(forces(3,:,16),forces(3,:,17)./forces(3,:,14),'d');
% set(gca, 'XDir','reverse')
% xlabel('Sinkage [m]')
% ylabel('Drawbar Pull [N] / Torque [N-m]')
% guide{4} = '\Delta b, W = 60 kg';
% guide{5} = '\Delta b, W = 90 kg';
% guide{6} = '\Delta b, W = 120 kg';
% legend(guide);

% figure()
% hold on
% for o = 1:length(slip_range)
% plot(forces(o,:,16),forces(o,:,14))
% end
% legend(num2str(slip_range'))
% xlabel('Sinkage [m]')
% ylabel('Torque [N-m]')
% set(gcf,'Position',[40 560 560 420])

figure()
hold on
for o = 1:length(Rw_range)
plot(slip_range(:),forces(:,o,15))
end
xlabel('Slip Ratio')
ylabel('Tractive Force [N]')
legend(num2str(Rw_range'));
set(gcf,'Position',[40 560 560 420])

figure()
hold on
for o = 1:length(Rw_range)
plot(slip_range(:),forces(:,o,17))
end
xlabel('Slip Ratio')
ylabel('Drawbar Pull [N]')
legend(num2str(Rw_range'));
set(gcf,'Position',[675 560 560 420])
xlim([0,1])
ylim([-100,100])

figure()
hold on
for o = 1:length(Rw_range)
plot(slip_range(:),forces(:,o,14))
end
xlabel('Slip Ratio')
ylabel('Torque [Nm]')
legend(num2str(Rw_range'));
set(gcf,'Position',[1300 560 560 420])
xlim([0,1])
ylim([0,100])

figure()
hold on
for o = 1:length(Rw_range)
plot(slip_range(:),forces(:,o,16))
end
xlabel('Slip Ratio')
ylabel('Sinkage [m]')
legend(num2str(Rw_range'));
set(gcf,'Position',[40 40 560 420])

figure()
hold on
for o = 1:length(Rw_range)
plot(slip_range(:),forces(:,o,11))
end
xlabel('Slip Ratio')
ylabel('Terrain Resistance [N]')
legend(num2str(Rw_range'));
set(gcf,'Position',[675 40 560 420])

