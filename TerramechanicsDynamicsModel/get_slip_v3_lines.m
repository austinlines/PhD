function [alpha, sigma] = get_slip_v3_lines(Romega, vx, vy, T)
%Function get_slip

%%% Changelog
% May 12, 2006
% Removed absolute values from several if test-cases -- Rw could be in
% theory larger in the negative direction than vx is in the positive
% direction, causing a false value.  (removals marked in code) -- JPM
% Added error messages for some cases for testing; commented out -- JPM
% Commented out the small_diff section (see comments) -- JPM

small_vx = 0.001;       % m/s
small_vy = 0.001;       % m/s
small_diff = 0.001;     % a difference of 1 cm/s is considered negligible here
if (abs(vx) < small_vx) || (abs(vy) < small_vy)
    alpha = 0;
else
    alpha = -atan2(vy, abs(vx)); % radians
end

v = sqrt(vx^2 + vy^2);
if (abs(vx) < small_vx)&&(abs(Romega) < small_vx)
    slip = 0;
% This small_diff case *seems* like it ought to be a good idea, but causes lots of problems in practice.    
% elseif (abs(Romega - vx) < small_diff)
%     slip = 0;
elseif abs(T) > 0           % Is a torque being applied?
    if (T > 0) && (vx > 0)
        if Romega > vx      % removed absolute values - jpm 5/12
            slip = 100*(Romega - vx)/(Romega);      % positive slip, traction
        elseif Romega < vx  
%             error('Romega < vx when T > 0')
%             slip = 100*(Romega - vx)/(vx);      % alines 2019-12-17, added
            slip = 0;                                     % alines 2019-12-17, commented out
        else
            slip = 0;       % This is the condition when torque is too low and T - w*R < 0
        end
    elseif (T < 0) && (vx > 0)
        slip = 100*(vx - Romega)/vx;
        slip = -abs(slip);              % positive slip braking
    elseif (T > 0) && (vx < 0)          % braking in reverse
        slip = 100*(vx - Romega)/vx;
        slip = abs(slip);
    elseif (T < 0) && (vx < 0)
        if Romega < vx % removed absolute values - jpm 5/12
            slip = 100*(Romega - vx)/(Romega);      % traction in reverse
        elseif Romega > vx
%             error('Romega > vx when T < 0');
            slip = 0;
        else
            slip = 0;                                 % This is the condition when torque is too low and T - w*R < 0
        end
        slip = -abs(slip);                              % traction in reverse
    % I'm not sure these next two are true.  The force drives down omega
    elseif (T > 0) && (abs(vx) < small_vx)
        slip = 100;
    elseif (T < 0) && (abs(vx) < small_vx)
        slip = -100;
    end
else   % no applied torque         
    if abs(vx) < eps  % vehicle not moving forward/backward
        if Romega > 0
            slip = 100;
        elseif Romega < 0
            slip = -100;
        end
    elseif vx > 0
        if Romega > vx
            slip = 100*(Romega - vx)/(Romega);  % kind of like traction
        elseif Romega < vx
            slip = -100*(vx - Romega)/(vx);  % braking
        else
            slip = 0;
        end
    elseif vx < 0
        if Romega < vx
            slip = -100*abs( (Romega - vx)/(Romega) );
        elseif Romega > vx
            slip = 100*abs( (vx - Romega)/(vx) );
        else
            slip = 0;
        end
    end
end
if slip > 100
    slip = 100;
elseif slip < -100
    slip = -100;
end
sigma = slip;
% disp(['Rw = ' num2str(Romega) ' vx = ' num2str(vx) ' T = ' num2str(T) ' sigma = ' num2str(sigma)])