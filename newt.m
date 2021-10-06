function [epsilon] = newt(gamma_meas,trials,theta_i,alpha)
%% Newton's method
%   gamma_meas is how we compute first epsilon/guess
% we assume 0 degree incidence ==> theta_i = 0
%  theta_t = asin(sin0/sqrt(epsilon)) = asin(0) = 0
% ==> gamma_meas(1+sqrt(epsilon)) = 1-sqrt(epsilon)
% ==> gamma_meas + gamma_meas*sqrt(epsilon) = 1 - sqrt(epsilon)
% ==> sqrt(epsilon)(gamma_meas+1) = 1 - gamma_meas
x1 = ((1-gamma_meas)/(1+gamma_meas))^2; % initial guess
x2 = x1*1.05;
for trial=1:trials
    theta_t1 = asin(sin(theta_i)/(sqrt(x1)));
    theta_t2 = asin(sin(theta_i)/(sqrt(x2)));
    y1 = gamma_meas - (cos(theta_t1)-cos(theta_i)*sqrt(x1))/(cos(theta_t1)+cos(theta_i)*sqrt(x1));
    y2 = gamma_meas - (cos(theta_t2)-cos(theta_i)*sqrt(x2))/(cos(theta_t2)+cos(theta_i)*sqrt(x2));
    if (abs(y1)<alpha)
        epsilon = x1;
        break
    elseif (abs(y2)<alpha)
        epsilon = x2;
        break
    else
        x3 = abs(((x2-x1)*(-y1)/(y2-y1))+x1);
        x1 = x3; x2 = x3*1.05;
    end
end
end

