function [xdot,acc,slip_terms] = f_of_x8_lines2(acc,x,u,ww,h,slip_terms,force_terms)
% Propagates state for 8 DOF nonlinear system: 

global g Lf Lr hcg R Rw mg m muf mur ms LfpLr Iw Iwf Iwr e track_f tr...
    ho Ixxs Ixzs Ixz Izz mhur mhuf mshcg mshcgg mshcgpho hfLrp hrLfm Fz
global Fzf_1 Fzr_1 tfo2 tro2 maxForce Long Lat Bw

% x = [vx vy r wfl wfr wrl wrr]

pos = @(x)x>0;
neg = @(x)x<0;
zer = @(x)x==0;

% State
vx = x(1);      % x velocity
vy = x(2);      % y velocity
r = x(3);       % yaw rate
wfl = x(4);     % angular velocity fl
wfr = x(5);     % angular velocity fr
wrl = x(6);     % angular velocity rl
wrr = x(7);     % angular velocity rr
ax = acc(1);    % x acceleration
ay = acc(2);    % y acceleration
rdot = acc(3);  % yaw acceleration

if vx == 0, vx = eps;  end

% Controls
Tfl = u(1);   Tfr = u(2);
Trl = u(3);   Trr = u(4);

% NORMAL FORCES AT EACH WHEEL
% Accelerations 
r2 = r^2;
ays  = ay + rdot*e;
axuf = ax - r2*Lf;
ayuf = ay + rdot*Lf;
axur = ax + r2*Lr;
ayur = ay - rdot*Lr;

% Longitudinal weight shift
Fzf2 = -(mhuf*axuf + mhur*axur)/LfpLr;

% Normal forces
Fzf = 0.5*(Fzf_1 + Fzf2);
Fzr = 0.5*(Fzr_1 - Fzf2);

% Calculate speeds at each tire
%   These calculations incorporate x and y velocity components and yaw (r)
%   times moment arm (either track_f, tr for front, rear track width or Lf, Lr for
%   front, rear distance from axle line to cg
vxfl = vx - r*track_f/2;
vyfl = vy + r*Lf;
vxfr = vx + r*track_f/2;
vyfr = vy + r*Lf;
vxrl = vx - r*tr/2;
vyrl = vy - r*Lr;
vxrr = vx + r*tr/2;
vyrr = vy - r*Lr;

% Retrieve slip and force values
alpha_fl = slip_terms(1);
alpha_fr = slip_terms(2);
alpha_rl = slip_terms(3);
alpha_rr = slip_terms(4);
sigma_fl = slip_terms(5);
sigma_fr = slip_terms(6);
sigma_rl = slip_terms(7);
sigma_rr = slip_terms(8);
Fxfl = force_terms(1);
Fxfr = force_terms(2);
Fxrl = force_terms(3);
Fxrr = force_terms(4);
Fyfl = force_terms(5);
Fyfr = force_terms(6);
Fyrl = force_terms(7);
Fyrr = force_terms(8);
Rxfl = force_terms(9);
Rxfr = force_terms(10);
Rxrl = force_terms(11);
Rxrr = force_terms(12);

% Normal forces
Tfl_r = force_terms(17);
Tfr_r = force_terms(18);
Trl_r = force_terms(19);
Trr_r = force_terms(20);

if length(force_terms) > 20
    Mzfl = force_terms(21);
    Mzfr = force_terms(22);
    Mzrl = force_terms(23);
    Mzrr = force_terms(24);
    Mz = Mzfl + Mzfr + Mzrl + Mzrr;
else
    Mzfl = 0;
    Mzfr = 0;
    Mzrl = 0;
    Mzrr = 0;
    Mz = 0;
end

%%
% Total lat./long. forces
Fx = Fxfl + Fxfr + Fxrr + Fxrl;
Fy = Fyfl + Fyfr + Fyrr + Fyrl;

%% Function evaluation
vxd = (Fx - Rxfl - Rxfr - Rxrl - Rxrr)/m;% + vy*r + ww(1) - 0.08*vx;   % 0.001
vyd = 0; %Fy/m - vx*r + ww(2);
Bt = 1;

% Yaw Rate
Mres = 100;
rd = (1/Izz)*((tr/2)*((Fxfr + Fxrr)-(Fxfl+Fxrl)) - Mres*r);
   
% Wheel Accelerations
wfld = ((Tfl - Bw*wfl - Tfl_r)/Iwf) + ww(4);    % kinetic friction
wfrd = ((Tfr - Bw*wfr - Tfr_r)/Iwf) + ww(5);    % kinetic friction
wrld = ((Trl - Bw*wrl - Trl_r)/Iwf) + ww(6);    % kinetic friction
wrrd = ((Trr - Bw*wrr - Trr_r)/Iwf) + ww(7);    % kinetic friction
 
acc = real([(vxd-vy*r) (vyd+vx*r) rd]);         % instantaneous accelerations
xdot(:,1) = real([vxd vyd rd wfld wfrd wrld wrrd]');
