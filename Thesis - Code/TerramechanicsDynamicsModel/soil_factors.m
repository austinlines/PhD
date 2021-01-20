function [Kc, Ky] = soil_factors(phi)

% Lookup table for bearing capacity factors of snow
phi_lookup = 0:2:50;
Nq_lookup = [1,1.22,1.49,1.81,2.21,2.69,3.29,4.02,4.92,6.04,7.44,9.19,11.40,14.21,...
    17.81,22.46,28.52,36.50,47.16,61.55,81.27,108.75,147.74,204.19,207.85,415.15];
Nc_lookup = [5.70,6.30,6.97,7.73,8.60,9.60,10.76,12.11,13.68,15.52,17.69,20.27,23.36,27.06,...
    31.61,37.16,44.04,52.64,63.53,77.50,95.66,119.67,151.95,196.22,258.29,347.51];
Ny_lookup = [0,0.2,0.4,0.6,0.9,1.2,1.7,2.3,3.0,3.9,4.9,5.8,7.8,11.7,...
    15.7,19.7,27.9,36.0,52.0,80.0,100.4,180.0,257.0,420.0,780.0,1153.2];

Nq = interp1(phi_lookup,Nq_lookup,rad2deg(phi));
Nc = interp1(phi_lookup,Nc_lookup,rad2deg(phi));
Ny = interp1(phi_lookup,Ny_lookup,rad2deg(phi));

% From Apostolopoulos, 2001
Kc = (Nc - tan(phi))*cos(phi)^2;
Ky = ((2*Ny)/tan(phi)+1)*cos(phi)^2;

end