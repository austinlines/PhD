function y  = shear_stress_lines(xt,theta1,~,kc,b,kphi,n,r,phi,c,K,Kr,Kw,slip,flag,model,km,zm,kp1,kp2,kz1,kz2,~,~,shear_model,z_offset)

% Function for normal stress integration from theta1 to thetam
z = r*(cos(xt) - cos(theta1)) + z_offset;
if model ==1
    sigma = ((kc/b + kphi)*(r*(cos(xt) - cos(theta1))).^n);
elseif model == 2
    
    sigma = zeros(size(z));
    sigma(z<zm) = km*zm*(-log(1-(z(z<zm)/zm)));
    x1 = zm-0.0002;
    x2 = zm-0.0001;
    y1 = km*zm*(-log(1-(x1/zm)));
    y2 = km*zm*(-log(1-(x2/zm)));
    slope = (y2-y1)/(x2-x1);
    if sum(z>=zm)>1
        sigma(z>=zm) = (z(z>=zm)-x1)*slope+y2;
    end
elseif model == 3
    sigma = (kp1 + kp2*b)*(-log(1-(z/(kz1+(kz2/b)))));
elseif model == 4
    sigma = ((kc + kphi*b)*((r/b)*(cos(xt) - cos(theta1))).^n);
end

j = r*(theta1 - xt - (1-slip).*(sin(theta1) - sin(xt)));

if shear_model == 2
    s = Kr*(1+(1/(Kr*(1-1/exp(1)))-1).*exp(1-j/Kw)).*(1-exp(-j/Kw));
else
    s = 1 - exp(-j/K);
end

tau = (c + sigma*tan(phi)).*s;

if flag == 'normal'
    y = tau.*sin(xt);
elseif flag == 'resist'
    y = tau.*cos(xt);
elseif flag == 'resis2'
    y = tau.*cos(xt);
elseif flag == 'torque'
    y = tau;
end


