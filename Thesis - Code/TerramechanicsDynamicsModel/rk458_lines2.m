function [yout,acc,slip_terms] = rk458_lines(acc,y,u,h,w,slip_terms,force_terms)
global beta gam

n=max(size(y));
if nargin < 6
    w = zeros(max(size(y)),1);
end

% Initialization
slipxdot = zeros(n,6);
% Compute the derivatives
[xdot(:,1),acc] = f_of_x8_lines2(acc,y,u,w,h,slip_terms,force_terms);
ypd = y + h*xdot(:,1)*0.25;
[xdot(:,2),acc] = f_of_x8_lines2(acc,ypd,u,w,h,slip_terms,force_terms);
ypd = y + h*xdot(:,1:2)*beta(1:2,2);
[xdot(:,3),acc] = f_of_x8_lines2(acc,ypd,u,w,h,slip_terms,force_terms);
ypd = y + h*xdot(:,1:3)*beta(1:3,3);
[xdot(:,4),acc] = f_of_x8_lines2(acc,ypd,u,w,h,slip_terms,force_terms);
ypd = y + h*xdot(:,1:4)*beta(1:4,4);
[xdot(:,5),acc] = f_of_x8_lines2(acc,ypd,u,w,h,slip_terms,force_terms);
ypd = y + h*xdot(:,1:5)*beta(1:5,5);
[xdot(:,6),acc] = f_of_x8_lines2(acc,ypd,u,w,h,slip_terms,force_terms);

% Estimate the error and the acceptable error
yout = y + h*xdot*gam;
[~,acc,slip_terms] = f_of_x8_lines2(acc,yout,u,w,h,slip_terms,force_terms);