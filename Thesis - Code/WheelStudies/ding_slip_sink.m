function [delta_z1] = ding_slip_sink(thetam,theta2,phi,slip)

global Rw

if slip < 0.01
    delta_z1 = 0;
else
    alpha0 = 2*phi - thetam;
    Xc = (pi()/4) - (phi/2);
    
    x_A = Rw*sin(thetam);
    y_A = Rw*(1-cos(thetam));
    
    x_E = Rw*sin(theta2);
    y_E = Rw*(1-cos(theta2));
    
    x_O = (y_A - y_E + x_A*tan(alpha0) + x_E*tan(Xc))/(tan(alpha0)+tan(Xc));
    y_O = y_A + tan(alpha0)*(x_A - x_O);
    
    rho0 = sqrt((x_A - x_O)^2 + (y_A - y_O)^2);
    
    syms alpha_A
    eqn = Rw*(1-slip)*thetam == x_O + rho0*exp((alpha_A - alpha0)*tan(phi))*cos(alpha_A);
    alpha_A = vpasolve(eqn,alpha_A);%,'Real',true);
    alphaA = double(alpha_A);
    
    delta_z1 = rho0*exp((alphaA - alpha0)*tan(phi))*sin(alphaA) - y_O;
    if delta_z1 > 0
        delta_z1 = 0;
    end
end

end
