function y  = normal_stress_lines(xt,theta1,~,kc,b,kphi,n,r,flag,model,km,zm,kp1,kp2,kz1,kz2,z_offset)

% Function for normal stress integration from theta1 to thetam
if (cos(xt) - cos(theta1)) < 0
    disp('neg');
end

z = r*(cos(xt) - cos(theta1))+z_offset;
if model == 1
    y = ((kc/b + kphi)*(r*(cos(xt) - cos(theta1))).^n);
elseif model == 2
    y = zeros(size(z));
    y(z<zm) = km*zm*(-log(1-(z(z<zm)/zm)));
    x1 = zm-0.0002;
    x2 = zm-0.0001;
    y1 = km*zm*(-log(1-(x1/zm)));
    y2 = km*zm*(-log(1-(x2/zm)));
    slope = (y2-y1)/(x2-x1);
    if sum(z>=zm)>1
        y(z>=zm) = (z(z>=zm)-x1)*slope+y2;
    end
elseif model == 3
    y = (kp1 + kp2*b)*(-log(1-(z/(kz1+(kz2/b)))));
elseif model == 4
    y = ((kc + kphi*b)*((r/b)*(cos(xt) - cos(theta1))).^n);
end


if flag == 'normal'
    y = y.*cos(xt);
elseif flag == 'resist'
    y = y.*sin(xt);
elseif flag == 'torque'
end    
    


