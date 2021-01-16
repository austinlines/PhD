clear
close all
load('p_mult.mat')
q = 1;

[filename,pathname] = uigetfile('.mat');
if filename==0
  return
end
load([pathname filename]);
start2 = start;
stop2 = stop;

sustain = 3;

% Number of Hypotheses
J_total = size(p,2);

while q < 5

% System model
% x_k = Phi*x_k-1 + Gamma*u_k-1 + Lambda*w_k-1
% z_k = H*x_k + n

tic
H = 1;                      % Measurement to Y translation
P = ones(J_total,1);        % Process Covariance = E([x-x^][x-x^]^T)
Lmb = .1;                   % Disturbance to X translation
Q = 20;                     % Disturbance Covariance
R = 1;                      % Measurement Error Covariance


%%
start = 1;

if q == 1
    u(1,:) = e.T_fl(start:end);
    u(2,:) = e.wfl(start:end);
    z = e.v_X(start:end);
elseif q == 2
    u(1,:) = e.T_fr(start:end);
    u(2,:) = e.wfr(start:end);
    z = e.v_X(start:end);
elseif q == 3
    u(1,:) = e.T_rl(start:end);
    u(2,:) = e.wrl(start:end);
    z = e.v_X(start:end);
elseif q == 4
    u(1,:) = e.T_rr(start:end);
    u(2,:) = e.wrr(start:end);
    z = e.v_X(start:end);
end

nt = size(u,2)-1;

minimum = 1/(1.1*J_total);

%% Initialize Filter
xhat = zeros(J_total,1);
posterior = ones(J_total,1)/J_total;

for k = 1:nt
    posterior(posterior<minimum) = minimum;
    
    min_count = sum(posterior==minimum);
    remainder = 1 - min_count*minimum;
    orig_rem = sum(posterior(posterior~=minimum));
    for i = 1:J_total
        if posterior(i) == minimum
            prior(i,1) = minimum;
        else
            prior(i,1) = remainder*(posterior(i)/orig_rem);
        end
    end
    
    % Resets prior to be a uniform distribution each time
    prior = ones(J_total,1)/J_total;

    for j = 1:size(p,2)
        % System model
        Phi = p(1,j);
        Gamma = p(2:3,j)';
        
        % Kalman filter
        xhatm = Phi*xhat(j) + Gamma*u(:,k);        % State estimate extrapolation - no u
        Pm = Phi*P(j)*Phi' + Lmb*Q*Lmb';           % Covariance extrapolation
        K = Pm*H'*inv(H*Pm*H' + R);                % filter gain
        xhat(j) = xhatm + K*(z(k) - H*xhatm);      % update state
        P(j) = inv(inv(Pm) + H'*inv(R)*H);         % update covariance
        
        % Measurement residual and residual covariance matrix
        r(j) = z(k) - H*xhatm;
        S(j) = H*P(j)*H + R;
         
        likelihood(j) =(1/sqrt(2*pi*S(j)))*exp(-(r(j)^2)/(2*S(j)));
        numerator(j) = likelihood(j)*prior(j);

        % Store results
        yhat(k,j) = xhat(j);

    end
    if sum(numerator) == 0
        posterior = zeros(size(posterior));
    else
        posterior = numerator/sum(numerator);
    end
    p_est(:,k) = sum(posterior.*p,2);
    index(k) = min(find(posterior == max(posterior)));
    post_sum1(k) = sum(posterior(1,1:J_nom))/J_nom;
    post_sum2(k) = sum(posterior(1,(J_nom+1):J_total))/(J_total-J_nom);
    xhat_exp(k) = sum(posterior.*xhat');
end
toc

%%
if q == 1
    index_fl = index;
    post_sum1_fl = post_sum1;
    post_sum2_fl = post_sum2;
elseif q == 2
    index_fr = index;
    post_sum1_fr = post_sum1;
    post_sum2_fr = post_sum2;
elseif q == 3
    index_rl = index;
    post_sum1_rl = post_sum1;
    post_sum2_rl = post_sum2;
elseif q == 4
    index_rr = index;
    post_sum1_rr = post_sum1;
    post_sum2_rr = post_sum2;
end

q = q+1;
end

%%
stop_all = length(e.t);

fig = figure();
subplot(5,1,1)
plot(e.t(1:stop_all),e.v_X(1:stop_all),'k');
ylabel({'Longitudinal';'Velocity [m/s]'})
xline(start2/10-.1);
xline(stop2/10-.1);
xline(detect_imm/10-.1,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2/10-.1,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3/10-.1,'--r');
end

subplot(5,1,2)
title('Front Left')
h1 = plot(e.t(1:stop_all),e.wfl(1:stop_all));
ylim([0,10])
yyaxis right
hold on
h2 = plot(e.t(1:stop_all),e.T_fl(1:stop_all),'-.');
h3 = plot(e.t(1:stop_all),e.slip_fl(1:stop_all),':','LineWidth',1.5,'Color',[0.4 0.8 0.3]);
ylabel(' ')
ylim([0,100])
xline(start2/10-.1);
xline(stop2/10-.1);
xline(detect_imm/10-.1,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2/10-.1,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3/10-.1,'--r');
end
legend('Location','northwest')
lgd_fr = legend([h1,h2,h3],{'wheel speed' 'torque' 'slip'});
title(lgd_fr,'Front Left')

subplot(5,1,3)
h1 = plot(e.t(1:stop_all),e.wfr(1:stop_all));
ylim([0,10])
yyaxis right
hold on
h2 = plot(e.t(1:stop_all),e.T_fr(1:stop_all),'-.');
h3 = plot(e.t(1:stop_all),e.slip_fr(1:stop_all),':','LineWidth',1.5,'Color',[0.4 0.8 0.3]);
ylim([0,100])
xline(start2/10-.1);
xline(stop2/10-.1);
xline(detect_imm/10-.1,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2/10-.1,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3/10-.1,'--r');
end
legend('Location','northwest')
lgd_fr = legend([h1,h2,h3],{'wheel speed' 'torque' 'slip'});
title(lgd_fr,'Front Right')

subplot(5,1,4)
h1 = plot(e.t(1:stop_all),e.wrl(1:stop_all));
ylim([0,10])
yyaxis right
hold on
h2 = plot(e.t(1:stop_all),e.T_rl(1:stop_all),'-.');
h3 = plot(e.t(1:stop_all),e.slip_rl(1:stop_all),':','LineWidth',1.5,'Color',[0.4 0.8 0.3]);
ylim([0,100])
xline(start2/10-.1);
xline(stop2/10-.1);
xline(detect_imm/10-.1,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2/10-.1,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3/10-.1,'--r');
end
legend('Location','northwest')
lgd_fr = legend([h1,h2,h3],{'wheel speed' 'torque' 'slip'});
title(lgd_fr,'Rear Left')

subplot(5,1,5)
title('rear right')
h1 = plot(e.t(1:stop_all),e.wrr(1:stop_all));
ylim([0,10])
yyaxis right
hold on
h2 = plot(e.t(1:stop_all),e.T_rr(1:stop_all),'-.');
h3 = plot(e.t(1:stop_all),e.slip_rr(1:stop_all),':','LineWidth',1.5,'Color',[0.4 0.8 0.3]);
ylim([0,100])
xline(start2/10-.1);
xline(stop2/10-.1);
xline(detect_imm/10-.1,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2/10-.1,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3/10-.1,'--r');
end
legend('Location','northwest')
lgd_fr = legend([h1,h2,h3],{'wheel speed' 'torque' 'slip'});
title(lgd_fr,'Rear Right')

set(gcf,'Position',[1120 65 800 900])
suplabel('Time [s]','x',[.25 .11 .55 .85]);
suplabel({'Wheel Speed';'[rad/s]'},'y',[.14 0 .55 .85]);
suplabel({'Torque [Nm]','Slip [%]'},'yy',[.35 0 .55 .85]);

stop_dist = e.dist(stop2)-e.dist(start2);
stop_time = e.t(stop2)-e.t(start2);
detect_dist = e.dist(detect_imm) - e.dist(start2);
detect_time = e.t(detect_imm) - e.t(start2);
dist_prop = (detect_dist/stop_dist)*100;
t_prop = (detect_time/stop_time)*100;

disp(['Stop Distance = ',num2str(stop_dist),' meters']);
disp(['Detect Distance = ',num2str(detect_dist),' meters']);
disp(['Detect Dist % = ',num2str(dist_prop),' %']);
disp(['Stop Time = ',num2str(stop_time),' sec']);
disp(['Detect Time = ',num2str(detect_time),' sec']);
disp(['Detect Time % = ',num2str(t_prop),' %']);
disp(['V at detect = ',num2str(e.v_X(detect_imm))]);
if J_total == 20
    disp('multi2')
elseif J_total == 27
    disp('multi3')
end


%%
sustain2 = 1;
diag1 = (post_sum1_fl + post_sum1_rr - post_sum2_fl - post_sum2_rr);
diag2 = (post_sum1_fr + post_sum1_rl - post_sum2_fr - post_sum2_rl);
front1 = (post_sum1_fl + post_sum1_fr - post_sum2_fl - post_sum2_fr);
rear1 = (post_sum1_rl + post_sum1_rr - post_sum2_rl - post_sum2_rr);

figure()
hold on
h1 = plot(movmean(diag1,[sustain2,0]));
h2 = plot(movmean(diag2,[sustain2,0]));
h3 = plot(movmean(front1,[sustain2,0]));
h4 = plot(movmean(rear1,[sustain2,0]));
legend('fl - rr','fr - rl','front','rear');
xline(start2);
xline(stop2);
xline(detect_imm,'--r');
if exist('detect_imm2')==1
    xline(detect_imm2,'--r');
end
if exist('detect_imm3')==1
    xline(detect_imm3,'--r');
end
grid on

